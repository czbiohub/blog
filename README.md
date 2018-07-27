# CZ Biohub DS Blog

This is the [Jekyll](https://github.com/jekyll/jekyll) project for the [CZ Biohub](https://www.czbiohub.org) Data Sciences group blog. It's based on the [Jekyll Now](http://www.jekyllnow.com/) project created by Barry Clark, with a couple bug-fixes and some add-ons:

- Posts can use [MathJax](https://www.mathjax.org/) by putting `mathjax: true` in the header section.
- Blocks of code can be hidden (showable via a button) with `codehide: true`.

For any other info about the framework, check out the Jekyll Now site and [repository](https://github.com/barryclark/jekyll-now).

### Making a post from a Jupyter notebook

Jupyter's `nbconvert` tool can convert notebooks into nicely-formatted Markdown. This is almost exactly what we need for a blogpost, except that a) it doesn't include the metadata header that Jekyll expects, and b) it doesn't put images in the right location relative to where the post will be.

To remedy those two issues, this repository includes a simple `jekyll.tpl` template which extends `nbconvert`'s `markdown.tpl` and makes the needed tweaks. Once your notebook is ready to be converted to Markdown, run the following:

```
jupyter nbconvert --to markdown --template path/to/this/repo/jekyll.tpl  YYYY-MM-DD-Your-Post-Title.ipynb
```

Note that your post filename needs to have the correct date format, or Jekyll will not publish it. After converting, move the resulting `.md` file to `_posts` and the folder of images to `images`. Make sure that things look right by running `jekyll serve --future`. If the future flag is not set, posts with upcoming dates will not be displayed on the site.
