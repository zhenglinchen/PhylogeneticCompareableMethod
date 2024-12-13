MoransI.R and ESRRV.R are the scipts(functions) to select the reasonable phylogenetic eigenvectors according mMorI and ESRRV method.



mMorI method: it is a stepwise algorithm, Firstly, it indicates the interesting trait as dependent varaible, and adds new phylogenetic eigenvectors as independent variables to the previous model, calculates the residuals' moran's I value, and then select the eigenvector that minimize the residuals' moran's I value in this round and keep it in the regression model, redo this process untile the residuals' moran's I value decreased to zero (or show no significance to zero).

ESRRV method: select the eigenvectors that significantly correlated to dependent variable.
