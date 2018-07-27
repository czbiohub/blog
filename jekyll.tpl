{% extends 'markdown.tpl' %}

{% block header -%}
---
layout: post
mathjax: true
codehide: true
title: {{ (resources['metadata']['name'][11:] or resources['metadata']['name']) | replace("-", " ") }}
---
{% endblock header %}

{% block data_svg -%}
![svg](/images/{{ output.metadata.filenames['image/svg+xml'] | path2url }})
{%- endblock data_svg %}

{% block data_png -%}
![png](/images/{{ output.metadata.filenames['image/png'] | path2url }})
{%- endblock data_png %}

{% block data_jpg -%}
![jpeg](/images/{{ output.metadata.filenames['image/jpeg'] | path2url }})
{%- endblock data_jpg %}
