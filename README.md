# E-vident: elucidating sampling effort for microbial analysis studies

A daunting question that confronts microbial ecologists when designing their studies is how many sequences and samples are sufficient to observe a specific effect size. An effect size is defined as the difference between groups of samples, such as the separation in β diversity between samples within groups of a study. Here we present E-vident, an interactive web application that allows us to graphically and statistically determine how many samples and sequences would have been necessary in order to find the effects reported in published studies. Results are presented in jackknifed plots, with multiple iterations to show statistical differences in β diversity. Using this tool, researchers can focus the sampling effort of new studies by relying on previously published data to capture specific differences between samples using statistical certainty.

Evident is a web-based software tool with an interactive user interface, implemented in HTML, Web Graphics Library, mod_python and Quantitative Insights Into Microbial Ecology ([QIIME][]). The interface of Evident (Fig. 1) is comprised of: (i) the selction of parameters (i.e., study of interest, number of sequences per sample, the number of samples to use, and the number of iterations); (ii) the kind of visualizations to generate (Demo PCoA shows the original results from the study, and PCoA recalculates the study using the user-defined parameters in (i); (iii) the WebGL plot display.

![Imgur](http://i.imgur.com/seMQ0.png)

*Figure 1*. E-vident GUI. Selectors for study, sequences per sample, number of subjects, samples per subject and number of iteration; 2. Analysis Menu: Demo (view original study results), PCoA plots and alpha rarefaction plots; C) Output area showing the webGL PCoA plot. B. Effect of rarefaction on β diversity; B1.

[QIIME]: https://github.com/qiime/qiime
