Reorganized in the style of a python package.

The folder that's currently "Analysis" should be eventually renamed whatever we want the top-level package to be called.

ViualStats will eventually be moved to bin. However, I know it is getting some use at the moment, and I don't want to mess up paths for anyone. That is, the folder "Rendering" is in the wrong place at the moment. And, Vstats isn't OOP but will be soon. It will stop breaking the pattern once it's rewritten as an OOP package.

The folder "bin" contains all and only scripts that can be used as bash utilities. Eventually, all of these scripts will reference the source contained in Analysis.

"Notebooks" contains iPython notebooks that will eventually all use the source in Analysis.
