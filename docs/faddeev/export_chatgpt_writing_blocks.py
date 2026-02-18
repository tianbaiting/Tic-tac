#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 功能：从 chatgpt.html 中提取写入块（data-writing-block=true）内容，优先写入块，备用 code block（尤其是 markdown 代码块）。每个块保存为一个 .md 文件，文件名包含索引、块 ID 和标题（如果有）。如果没有写入块，则根据 code block 的语言或内容特征导出 markdown 代码块。使用 BeautifulSoup 解析 HTML，支持命令行参数指定输入文件、输出目录和数据源类型。 
import argparse
import re
from pathlib import Path

try:
    from bs4 import BeautifulSoup, NavigableString, Tag
except ImportError as e:
    raise SystemExit(
        "需要 BeautifulSoup：pip install beautifulsoup4\n"
        f"ImportError: {e}"
    )

MATH_BLOCK_RE = re.compile(r"^\[\s*\n(.*?)\n\s*\]$", re.DOTALL)
MARKDOWN_LANGS = {"markdown", "md", "mdx", "mkd", "mkdn", "mdown"}
COPY_BUTTON_LABELS = {"copy", "copy code", "复制", "复制代码"}
LANGUAGE_RE = re.compile(r"language-([A-Za-z0-9_+.-]+)")
SHORT_LANG_RE = re.compile(r"^[A-Za-z0-9_+.-]{1,24}$")


def normalize_ws(s: str) -> str:
    s = s.replace("\r\n", "\n").replace("\r", "\n")
    return s.replace("\u00a0", " ")


def safe_filename(s: str) -> str:
    s = s.strip()
    s = re.sub(r"[\\/:*?\"<>|]+", "_", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s[:80] if s else "untitled"


def get_text(node: Tag) -> str:
    return normalize_ws(node.get_text(separator="", strip=False))


def escape_md_text(s: str) -> str:
    return s


def md_from_inline(node) -> str:
    if isinstance(node, NavigableString):
        return escape_md_text(str(node))

    if not isinstance(node, Tag):
        return ""

    name = node.name.lower()

    if name == "br":
        return "\n"

    if name in ("span",):
        return "".join(md_from_inline(c) for c in node.children)

    if name in ("strong", "b"):
        inner = "".join(md_from_inline(c) for c in node.children).strip()
        return f"**{inner}**" if inner else ""

    if name in ("em", "i"):
        inner = "".join(md_from_inline(c) for c in node.children).strip()
        return f"*{inner}*" if inner else ""

    if name == "code":
        inner = get_text(node).replace("`", "\\`")
        return f"`{inner}`"

    if name == "a":
        href = node.get("href", "").strip()
        text = "".join(md_from_inline(c) for c in node.children).strip() or href
        return f"[{text}]({href})" if href else text

    return "".join(md_from_inline(c) for c in node.children)


def md_from_list_item(node: Tag, indent: int, ordered: bool, number: int) -> str:
    bullet = f"{number}. " if ordered else "- "
    prefix = " " * indent + bullet

    main_chunks = []
    sublists = []
    for c in node.children:
        if isinstance(c, Tag) and c.name.lower() in ("ul", "ol"):
            sublists.append(c)
        else:
            main_chunks.append(c)

    main_text = "".join(md_from_inline(c) for c in main_chunks).strip()
    if main_text:
        line = prefix + main_text.replace("\n", "\n" + " " * (indent + 2))
        out = line + "\n"
    else:
        out = prefix.rstrip() + "\n"

    for sub in sublists:
        out += md_from_block(sub, indent=indent + 2)

    return out


def md_from_block(node, indent: int = 0) -> str:
    if isinstance(node, NavigableString):
        s = str(node)
        return "" if s.strip() == "" else escape_md_text(s)

    if not isinstance(node, Tag):
        return ""

    name = node.name.lower()

    if name in ("div", "section"):
        if node.get("contenteditable") == "false" and node.find("hr"):
            return "\n\n---\n\n"
        parts = []
        for c in node.children:
            part = md_from_block(c, indent=indent)
            if part:
                parts.append(part)
        return "".join(parts)

    if name in ("h1", "h2", "h3", "h4", "h5", "h6"):
        level = int(name[1])
        text = "".join(md_from_inline(c) for c in node.children).strip()
        return f"\n{'#' * level} {text}\n\n" if text else ""

    if name == "hr":
        return "\n\n---\n\n"

    if name == "p":
        text = "".join(md_from_inline(c) for c in node.children)
        text = normalize_ws(text).strip("\n")
        if text.strip() == "":
            return "\n"
        m = MATH_BLOCK_RE.match(text.strip())
        if m:
            inner = m.group(1).strip("\n")
            return f"\n$$\n{inner}\n$$\n\n"
        return f"\n{text}\n\n"

    if name == "blockquote":
        inner = "".join(md_from_block(c, indent=indent) for c in node.children)
        inner = normalize_ws(inner).strip()
        if not inner:
            return ""
        lines = [ln for ln in inner.split("\n") if ln.strip()]
        quoted = "\n".join("> " + ln for ln in lines)
        return f"\n{quoted}\n\n"

    if name in ("ul", "ol"):
        parts = []
        items = node.find_all("li", recursive=False)
        ordered = name == "ol"
        for i, li in enumerate(items, start=1):
            parts.append(md_from_list_item(li, indent=indent, ordered=ordered, number=i))
        return "".join(parts) + ("\n" if parts else "")

    if name == "li":
        return md_from_list_item(node, indent=indent, ordered=False, number=1)

    if name == "pre":
        code = node.find("code")
        if code:
            code_text = get_text(code).strip("\n")
            lang = ""
            cls = " ".join(code.get("class", []))
            m = re.search(r"language-([A-Za-z0-9_+-]+)", cls)
            if m:
                lang = m.group(1)
            return f"\n```{lang}\n{code_text}\n```\n\n"
        pre_text = get_text(node).strip("\n")
        return f"\n```\n{pre_text}\n```\n\n"

    return "".join(md_from_block(c, indent=indent) for c in node.children)


def extract_writing_blocks(html: str):
    soup = BeautifulSoup(html, "html.parser")
    blocks = soup.select('div[data-writing-block="true"]')
    results = []

    for idx, wb in enumerate(blocks, start=1):
        wb_id = wb.get("id", f"writing-block-{idx:03d}")
        pm = wb.select_one("div.ProseMirror") or wb.select_one("div.writing-block-editor")
        if not pm:
            continue

        title_tag = pm.find(["h1", "h2", "h3"])
        title = get_text(title_tag).strip() if title_tag else ""

        md = md_from_block(pm, indent=0).strip()
        if not md:
            continue

        results.append(
            {
                "index": idx,
                "id": wb_id,
                "title": title,
                "markdown": md + "\n",
            }
        )

    return results


def infer_language_from_header(pre: Tag) -> str:
    for raw in pre.stripped_strings:
        label = raw.strip()
        if not label:
            continue
        lower = label.lower()
        if lower in COPY_BUTTON_LABELS:
            continue
        if not SHORT_LANG_RE.fullmatch(label):
            continue
        return lower
    return ""


def infer_language(pre: Tag, code: Tag) -> str:
    for node in (code, pre):
        class_text = " ".join(node.get("class", []))
        m = LANGUAGE_RE.search(class_text)
        if m:
            return m.group(1).lower()
        for attr in ("data-language", "data-lang", "lang"):
            value = str(node.get(attr, "")).strip().lower()
            if value:
                return value
    return infer_language_from_header(pre)


def looks_like_markdown(text: str) -> bool:
    lines = [line.rstrip() for line in normalize_ws(text).split("\n")]
    non_empty = [line for line in lines if line.strip()]
    if len(non_empty) < 2:
        return False

    checks = [
        any(re.match(r"^#{1,6}\s+\S", line) for line in non_empty),
        any(re.match(r"^\s*([-*+]\s+\S|\d+\.\s+\S)", line) for line in non_empty),
        any(re.match(r"^>\s+\S", line) for line in non_empty),
        any(re.match(r"^```", line) for line in non_empty),
        any(re.search(r"\[[^\]]+\]\([^)]+\)", line) for line in non_empty),
        any(re.match(r"^\|.+\|$", line) for line in non_empty),
    ]
    return sum(checks) >= 2


def extract_code_blocks(html: str, only_markdown: bool):
    soup = BeautifulSoup(html, "html.parser")
    results = []

    for source_index, pre in enumerate(soup.find_all("pre"), start=1):
        code = pre.find("code")
        if not code:
            continue

        content = normalize_ws(code.get_text(separator="", strip=False)).strip("\n")
        if not content.strip():
            continue

        language = infer_language(pre, code)
        is_markdown_lang = language in MARKDOWN_LANGS
        is_markdown_like = looks_like_markdown(content)

        if only_markdown and not (is_markdown_lang or is_markdown_like):
            continue

        results.append(
            {
                "index": len(results) + 1,
                "source_index": source_index,
                "language": language,
                "content": content + "\n",
            }
        )

    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="导出 chatgpt.html 内容到多个 .md 文件（优先写入块，备用 code block）。"
    )
    parser.add_argument("--input", default="chatgpt.html", help="输入 HTML 文件")
    parser.add_argument("--output-dir", default="md_blocks", help="输出目录")
    parser.add_argument(
        "--source",
        choices=("auto", "writing-blocks", "codeblocks"),
        default="auto",
        help="数据源：auto(默认优先写入块)、writing-blocks、codeblocks",
    )
    parser.add_argument(
        "--all-codeblocks",
        action="store_true",
        help="导出所有 code block（不只 markdown）",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        print(f"输入文件不存在：{in_path}")
        return

    html = in_path.read_text(encoding="utf-8")
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mode = args.source
    writing_items = []
    code_items = []
    used_code_fallback = False

    if mode in ("auto", "writing-blocks"):
        writing_items = extract_writing_blocks(html)
        if mode == "writing-blocks":
            if not writing_items:
                print("没有找到写入块 data-writing-block=true。")
                return
            for it in writing_items:
                fname = f"{it['index']:03d}_{safe_filename(it['id'])}"
                if it["title"]:
                    fname += f"__{safe_filename(it['title'])}"
                path = out_dir / f"{fname}.md"
                path.write_text(it["markdown"], encoding="utf-8")
            print(f"写出 {len(writing_items)} 个 Markdown 文件到：{out_dir.resolve()}")
            return

    if not writing_items and mode in ("auto", "codeblocks"):
        only_markdown = not args.all_codeblocks
        code_items = extract_code_blocks(html, only_markdown=only_markdown)
        if not code_items and only_markdown:
            code_items = extract_code_blocks(html, only_markdown=False)
            used_code_fallback = bool(code_items)

        if mode == "codeblocks":
            if not code_items:
                print("没有找到 code block。请确认导出 HTML 结构。")
                return
            for it in code_items:
                fname = f"{it['index']:03d}_codeblock_{it['source_index']:03d}"
                if it["language"]:
                    fname += f"__{safe_filename(it['language'])}"
                path = out_dir / f"{fname}.md"
                path.write_text(it["content"], encoding="utf-8")
            if used_code_fallback:
                print("没有检测到明确 markdown 代码块，已自动导出全部 code block。")
            print(f"写出 {len(code_items)} 个 Markdown 文件到：{out_dir.resolve()}")
            return

    if writing_items:
        for it in writing_items:
            fname = f"{it['index']:03d}_{safe_filename(it['id'])}"
            if it["title"]:
                fname += f"__{safe_filename(it['title'])}"
            path = out_dir / f"{fname}.md"
            path.write_text(it["markdown"], encoding="utf-8")
        print(f"写出 {len(writing_items)} 个写入块 Markdown 文件到：{out_dir.resolve()}")
        return

    if code_items:
        for it in code_items:
            fname = f"{it['index']:03d}_codeblock_{it['source_index']:03d}"
            if it["language"]:
                fname += f"__{safe_filename(it['language'])}"
            path = out_dir / f"{fname}.md"
            path.write_text(it["content"], encoding="utf-8")
        if used_code_fallback:
            print("没有检测到明确 markdown 代码块，已自动导出全部 code block。")
        print(f"写出 {len(code_items)} 个 code block Markdown 文件到：{out_dir.resolve()}")
        return

    print("没有找到可导出的写入块或 code block。")


if __name__ == "__main__":
    main()
