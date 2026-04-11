# tests/test_build_codebase_index.py
import ast
import pathlib
import sys
import textwrap

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "scripts"))
try:
    import build_codebase_index as bci
except ImportError:
    pytest.skip("scripts/build_codebase_index.py not available", allow_module_level=True)


class TestShouldSkipPath:
    def test_archive_component_skipped(self):
        assert bci.should_skip_path(pathlib.Path("model/_archive/foo.py"))

    def test_archive_as_substring_not_skipped(self):
        # "archive_old" does NOT contain "_archive" as a path component
        assert not bci.should_skip_path(pathlib.Path("model/archive_old/foo.py"))

    def test_pycache_skipped(self):
        assert bci.should_skip_path(pathlib.Path("model/__pycache__/foo.pyc"))

    def test_pyc_extension_skipped(self):
        assert bci.should_skip_path(pathlib.Path("model/foo.pyc"))

    def test_normal_path_not_skipped(self):
        assert not bci.should_skip_path(pathlib.Path("model/deconvolution/cuopt_impl.py"))


class TestFormatSignature:
    def _parse(self, src):
        return ast.parse(textwrap.dedent(src)).body[0]

    def test_simple_args(self):
        assert bci.format_signature(self._parse("def foo(a, b): pass")) == "foo(a, b)"

    def test_default_value(self):
        assert bci.format_signature(self._parse("def foo(a, b=1.0): pass")) == "foo(a, b=1.0)"

    def test_return_annotation(self):
        assert bci.format_signature(self._parse("def foo(a) -> int: pass")) == "foo(a) -> int"

    def test_self_skipped(self):
        assert bci.format_signature(self._parse("def foo(self, a): pass")) == "foo(a)"

    def test_cls_skipped(self):
        assert bci.format_signature(self._parse("def foo(cls, a): pass")) == "foo(a)"

    def test_star_args(self):
        assert bci.format_signature(self._parse("def foo(*args, **kwargs): pass")) == "foo(*args, **kwargs)"

    def test_type_annotations_excluded_from_args(self):
        # type hints on args should NOT appear, but return type SHOULD
        assert bci.format_signature(
            self._parse("def foo(a: int, b: str = 'x') -> bool: pass")
        ) == "foo(a, b='x') -> bool"

    def test_positional_only_args(self):
        node = self._parse("def foo(a, /, b): pass")
        assert bci.format_signature(node) == "foo(a, b)"


class TestWriteIfChanged:
    def test_writes_new_file(self, tmp_path):
        f = tmp_path / "out.md"
        changed = bci.write_if_changed(f, "header\nbody line")
        assert changed is True
        assert f.read_text() == "header\nbody line"

    def test_skips_when_body_unchanged(self, tmp_path):
        f = tmp_path / "out.md"
        f.write_text("old_header\nbody line")
        changed = bci.write_if_changed(f, "new_header\nbody line")
        assert changed is False
        assert f.read_text() == "old_header\nbody line"  # original preserved

    def test_writes_when_body_changed(self, tmp_path):
        f = tmp_path / "out.md"
        f.write_text("old_header\nold body")
        changed = bci.write_if_changed(f, "new_header\nnew body")
        assert changed is True
        assert f.read_text() == "new_header\nnew body"

    def test_creates_parent_dirs(self, tmp_path):
        f = tmp_path / "nested" / "dir" / "out.md"
        bci.write_if_changed(f, "header\nbody")
        assert f.exists()


class TestGetModuleSummary:
    def test_docstring_extracted(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text('"""This is the module."""\n\nimport os\n')
        assert bci.get_module_summary(f) == "This is the module."

    def test_multiline_docstring_first_line_only(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text('"""First line.\n\nMore detail.\n"""\n')
        assert bci.get_module_summary(f) == "First line."

    def test_comment_block_fallback(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("# Main comment\n# second line\nimport os\n")
        assert bci.get_module_summary(f) == "Main comment second line"

    def test_comment_block_stops_at_blank_line(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("# First block\n\n# Second block\nimport os\n")
        assert bci.get_module_summary(f) == "First block"

    def test_no_summary_returns_none(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("import os\n\ndef foo(): pass\n")
        assert bci.get_module_summary(f) is None

    def test_syntax_error_returns_none(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("def foo(\n")
        assert bci.get_module_summary(f) is None


class TestExtractFullDetail:
    def test_public_function_extracted(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text('def run_qp(adata, lam=1.0):\n    """Run QP solver."""\n    pass\n')
        info = bci.extract_full_detail(f)
        assert not info["error"]
        assert len(info["functions"]) == 1
        assert info["functions"][0]["sig"] == "run_qp(adata, lam=1.0)"
        assert info["functions"][0]["doc"] == "Run QP solver."

    def test_private_function_skipped(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("def _helper(): pass\ndef public(): pass\n")
        info = bci.extract_full_detail(f)
        assert len(info["functions"]) == 1
        assert info["functions"][0]["sig"] == "public()"

    def test_class_public_methods_only(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text(
            "class Solver:\n"
            "    \"\"\"A solver.\"\"\"\n"
            "    def solve(self, x):\n"
            "        \"\"\"Solve it.\"\"\"\n"
            "        pass\n"
            "    def _internal(self):\n"
            "        pass\n"
            "    def __init__(self):\n"
            "        pass\n"
        )
        info = bci.extract_full_detail(f)
        assert len(info["classes"]) == 1
        cls = info["classes"][0]
        assert cls["name"] == "Solver"
        assert cls["doc"] == "A solver."
        assert len(cls["methods"]) == 1
        assert cls["methods"][0]["sig"] == "solve(x)"
        assert cls["methods"][0]["doc"] == "Solve it."

    def test_staticmethod_prefixed(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text(
            "class Foo:\n"
            "    @staticmethod\n"
            "    def make(): pass\n"
            "    @classmethod\n"
            "    def from_dict(cls, d): pass\n"
        )
        info = bci.extract_full_detail(f)
        method_sigs = [m["sig"] for m in info["classes"][0]["methods"]]
        assert "[static] make()" in method_sigs
        assert "[cls] from_dict(d)" in method_sigs

    def test_syntax_error_returns_error_dict(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("def foo(\n")
        info = bci.extract_full_detail(f)
        assert info["error"] is True
        assert info["functions"] == []
        assert info["classes"] == []


class TestIntegration:
    """Run the generator against a tiny fake project tree."""

    def _make_project(self, tmp_path):
        # Package with __init__.py
        pkg = tmp_path / "mypkg" / "core"
        pkg.mkdir(parents=True)
        (pkg / "__init__.py").write_text("")
        (pkg / "solver.py").write_text(
            '"""Solve things."""\n\n'
            'def run(data, lam=1.0):\n'
            '    """Run the solver."""\n'
            '    pass\n\n'
            'def _internal(): pass\n'
        )

        # Script dir (file-level only)
        scripts = tmp_path / "scripts"
        scripts.mkdir()
        (scripts / "run_benchmark.py").write_text(
            "# Benchmark runner script\nimport os\n"
        )

        return tmp_path

    def test_generates_tier2_with_function(self, tmp_path):
        proj = self._make_project(tmp_path)
        full_roots = [proj / "mypkg"]
        file_roots = [proj / "scripts"]
        result = bci.build_tier2(full_roots, file_roots, proj)
        assert "solver.py" in result
        assert "run(data, lam=1.0)" in result
        assert "_internal" not in result   # private excluded

    def test_generates_tier1_table(self, tmp_path):
        proj = self._make_project(tmp_path)
        full_roots = [proj / "mypkg"]
        result = bci.build_tier1(full_roots, proj)
        assert "| Subpackage |" in result
        assert "solver.py" in result
        assert "run" in result

    def test_write_if_changed_prevents_rewrite(self, tmp_path):
        proj = self._make_project(tmp_path)
        out = tmp_path / "docs" / "index.md"
        full_roots = [proj / "mypkg"]
        file_roots = [proj / "scripts"]
        tier2 = "<!-- git:abc -->\n" + bci.build_tier2(full_roots, file_roots, proj)
        bci.write_if_changed(out, tier2)
        mtime1 = out.stat().st_mtime
        # Write again with different hash but same body
        tier2_same_body = "<!-- git:xyz -->\n" + bci.build_tier2(full_roots, file_roots, proj)
        bci.write_if_changed(out, tier2_same_body)
        assert out.stat().st_mtime == mtime1  # not rewritten


class TestInstallHook:
    def _make_git_dir(self, tmp_path):
        hooks_dir = tmp_path / ".git" / "hooks"
        hooks_dir.mkdir(parents=True)
        return tmp_path

    def test_creates_new_hook(self, tmp_path):
        proj = self._make_git_dir(tmp_path)
        bci.install_hook(proj)
        hook = proj / ".git" / "hooks" / "pre-commit"
        assert hook.exists()
        assert "build_codebase_index" in hook.read_text()
        assert oct(hook.stat().st_mode)[-3:] == "755"

    def test_idempotent_if_already_installed(self, tmp_path):
        proj = self._make_git_dir(tmp_path)
        bci.install_hook(proj)
        content_after_first = (proj / ".git" / "hooks" / "pre-commit").read_text()
        bci.install_hook(proj)
        content_after_second = (proj / ".git" / "hooks" / "pre-commit").read_text()
        assert content_after_first == content_after_second

    def test_appends_to_existing_hook_without_exit(self, tmp_path):
        proj = self._make_git_dir(tmp_path)
        hook = proj / ".git" / "hooks" / "pre-commit"
        hook.write_text("#!/usr/bin/env bash\necho hello\n")
        bci.install_hook(proj)
        content = hook.read_text()
        assert "echo hello" in content
        assert "build_codebase_index" in content

    def test_inserts_before_exit_in_existing_hook(self, tmp_path):
        proj = self._make_git_dir(tmp_path)
        hook = proj / ".git" / "hooks" / "pre-commit"
        hook.write_text("#!/usr/bin/env bash\nset -e\necho hello\nexit 0\n")
        bci.install_hook(proj)
        content = hook.read_text()
        lines = content.splitlines()
        index_call = next(i for i, l in enumerate(lines) if "build_codebase_index" in l)
        exit_line = next(i for i, l in enumerate(lines) if l.strip().startswith("exit"))
        assert index_call < exit_line  # generator call comes before exit

    def test_appends_when_exit_is_guarded_inside_if_block(self, tmp_path):
        # Regression: install_hook previously inserted before indented `exit 1`
        # inside an if block, so the generator only ran in the error path.
        # The fix: only match top-level (unindented) exit lines.
        proj = self._make_git_dir(tmp_path)
        hook = proj / ".git" / "hooks" / "pre-commit"
        hook.write_text(
            "#!/usr/bin/env bash\n"
            "status=0\n"
            "if [ \"$status\" -ne 0 ]; then\n"
            "    echo 'commit aborted'\n"
            "    exit 1\n"   # indented — must NOT be used as insertion point
            "fi\n"
        )
        bci.install_hook(proj)
        content = hook.read_text()
        lines = content.splitlines()
        # generator block must come AFTER the closing fi, not inside the if block
        fi_line = next(i for i, l in enumerate(lines) if l.strip() == "fi")
        index_call = next(i for i, l in enumerate(lines) if "build_codebase_index" in l)
        assert index_call > fi_line  # generator is after fi, not inside the if block
