// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec reference_matrix(const arma::mat& X,
                           const arma::vec& ind) {

  int cell_count = X.n_cols;

  vec reference_scores(cell_count);
  reference_scores.zeros();

  for (int i = 0; i < cell_count; i++){
    double sum = 0;

    for(int j = 0; j < ind.size(); j++){
      sum += X.at(ind[j], i);
    }
    reference_scores[i] = sum;
  }
  return reference_scores;
}
