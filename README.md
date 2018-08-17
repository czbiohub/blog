# CZ Biohub DS Blog

This is the [Jekyll](https://github.com/jekyll/jekyll) project for the [CZ Biohub](https://www.czbiohub.org) Data Sciences group blog. It's based on the [Jekyll Now](http://www.jekyllnow.com/) project created by Barry Clark, with a couple bug-fixes and some add-ons:

- Posts can use [MathJax](https://www.mathjax.org/) by putting `mathjax: true` in the header section.
- Blocks of code can be hidden (showable via a button) with `codehide: true`.

For any other info about the framework, check out the Jekyll Now site and [repository](https://github.com/barryclark/jekyll-now). In particular, for local development see [here](https://github.com/barryclark/jekyll-now#local-development).

### Making a post from a Jupyter notebook

Jupyter's `nbconvert` tool can convert notebooks into nicely-formatted Markdown. This is almost exactly what we need for a blogpost, except that a) it doesn't include the metadata header that Jekyll expects, and b) it doesn't put images in the right location relative to where the post will be.

To remedy those two issues, this repository includes a simple `jekyll.tpl` template which extends `nbconvert`'s `markdown.tpl` and makes the needed tweaks. Once your notebook is ready to be converted to Markdown, run the following:

```
jupyter nbconvert --to markdown --template path/to/this/repo/jekyll.tpl  Your-Post-Title.ipynb
```

You should add `date: YYYY-MM-DD` to the post's frontmatter. If the date is in the future, Jekyll won't publish it by default: use `jekyll serve --future` to view it locally. After converting, move the resulting `.md` file to `_blog` and the folder of images to `images`.

### Local testing

You can serve a local version of the blog by running `jekyll serve` in the root directory of this repo, after installing the relevant [dependencies](https://github.com/barryclark/jekyll-now#local-development). An easy workflow is to make a new branch for each post, write and preview it locally, and push it for review. Merge it into master when you're ready for primetime.
