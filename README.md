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

### Local install

Getting everything installed can be a little tricky. Here's what I did:

 1. Install Ruby using [homebrew](https://brew.sh) and `brew install ruby`. Installing using anaconda didn't work for me due to compiler issues that I don't want to debug. OS X includes Ruby already, but it's an older version and Jekyll needs something more recent.
      * Because OS X has Ruby already, homebrew will not put this version on your PATH. If you want you can set that up in your shell profile, but I decided to set up a blog conda environment (below).
 2. (optional) Set up a conda environment so that you can `conda activate [env-name]` when you want the blog to work.
     * create a file `~/miniconda3/envs/[env-name]/etc/conda/activate.d/env_vars.sh` that adds the relevant directories to your path, and another, `~/miniconda3/envs/[env-name]/etc/conda/deactivate.d/env_vars.sh` to remove it:
   
```shell
#!/bin/sh

export PATH="/usr/local/opt/ruby/bin:/usr/local/lib/ruby/gems/2.6.0/bin:$PATH"
```

```shell
#!/bin/sh

export PATH=${PATH#/usr/local/opt/ruby/bin:/usr/local/lib/ruby/gems/2.6.0/bin:}
```

 3. Install the needed [Ruby dependencies](https://github.com/barryclark/jekyll-now#local-development) with `gem install github-pages`


### Local testing

You can serve a local version of the blog by running `jekyll serve` in the root directory of this repo, after getting set up (above). An easy workflow is to make a new branch for each post, write and preview it locally, and push it for review. Merge it into master when you're ready for primetime.

