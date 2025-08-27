# bio_utils
Sometimes you end up re-inventing the wheel because an existing solution doesn't do *quite* what you want. This is my collection of wheels.

It contains utility functions that I've written over the years for recurring tasks in bioinformatics. Much functionality is duplicated in libraries like [scikit-bio](https://scikit.bio/) or [biotite](https://www.biotite-python.org/latest/), and for larger or more polished projects, I recommend their IO interfaces over the ones here. However, when I only need a simple parser or want to display a plot just the way I like it, I often find myself reaching for these functions since they're lightweight and lack the dependencies or complex class hierarchies of larger libraries.

Use the following commands to install `bio_utils` into your local environment:

```
git clone
cd bio_utils
pip install -e .
```
