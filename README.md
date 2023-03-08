# BayesMultiLogit

This is an R package developed by Kyle R. McEvoy and Jared D. Fisher to implement a number of methods for performing Bayesian multinomial logistic regression using data-augmentation methods. The data augmentation method underlying multilogit and multilogit_ESS is outlined in a paper by Jared D. Fisher and Kyle R. McEvoy titled "Bayesian Multinomial Logistic Regression for Numerous Categories", currently a work in progress.  

The code for all functions using the Holmes-Held methods were written by us following the pseudo-code template provided in the paper:
"Bayesian Auxiliary Variable Models for Binary and Multinomial Regression" by Chris C. Holmes and Leonhard Held, Bayesian Analysis (2006).

The code for the Polya-Gamma functions was copied from the public repository https://github.com/jgscott/helloPG by James G. Scott who was adapting code
from the package https://github.com/jwindle/BayesLogit. The BayesLogit repository includes code written by Jesse Windle, Nicholas Polson, and James G. Scott.

The multivariate normal density function was copied from work by Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen on the Rcpp Gallery
located at https://gallery.rcpp.org/articles/dmvnorm_arma/.

# License   

All of the work in this package falls under the Gnu General Public License Version 3 found at https://www.gnu.org/licenses/gpl-3.0.en.html.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


