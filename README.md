# Arellano-Bond-LASSO

Here, you can find the R codes to implement the simulation study and application described in the article 'Arellano-Bond LASSO Estimator for Dynamic Linear Panel Models' by Chernozhukov V, Fernández-Val I, Huang C and Wang W. (arXiv preprint arXiv: 2402.00584: https://arxiv.org/abs/2402.00584).

The dataset for the application on COVID-19 spread and school policy effects is initially provided by Chernozhukov V, Kasahara H, and Schrimpf P on the GitHub repository https://github.com/ubcecon/covid-schools. 
We aggregate the observations at the weekly level to avoid spurious serial correlation coming from the moving averages. Counties with missing values are dropped to obtain a balanced panel dataset with N=2510 counties and T=32 weeks.
