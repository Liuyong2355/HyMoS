# HyMoS-Baltamatica 参赛演示文稿

Beamer 源稿采用 16:10 页面，共 12 页。模板及 FIM、NMG 内联 TikZ 图沿用指定的 `beamer_v1.3`，算法图、表格及文字均可编辑。

## 文件

- `HyMoS-Baltamatica-slides.tex`：演示文稿源文件。
- `HyMoS-Baltamatica-slides.pdf`：用于展示的 PDF。
- `smile_styles.tex`：衬线字体、蓝色标题横线及页脚。
- `figures/`：使用流程中的实际插件截图。

## 编译

在本目录运行两遍 XeLaTeX：

```text
xelatex -interaction=nonstopmode -halt-on-error HyMoS-Baltamatica-slides.tex
xelatex -interaction=nonstopmode -halt-on-error HyMoS-Baltamatica-slides.tex
```

需要包含 `ctexbeamer`、TikZ 和 `listings` 的 TeX Live。参考文献在对应页面列出，无需运行 BibTeX。已在 TeX Live 2025 中编译。
