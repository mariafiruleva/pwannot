// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat cumulative_matrix(const arma::mat& X,
                            const arma::uvec& pathway_sizes,
                            const arma::mat& reference_scores,
                            int sampling_size = 100) {
  int genes_count = X.n_rows;
  int cell_count = X.n_cols;
  int max_size = pathway_sizes.max();

  arma::mat sample(pathway_sizes.size(), cell_count);
  arma::mat sum(max_size, cell_count);
  arma::mat results(pathway_sizes.size(), cell_count);
  results.zeros();
  for (int k = 0; k < sampling_size; k++) {
    uvec random_samples = randi<uvec>(max_size, distr_param(0, genes_count-1));
    sum = arma::cumsum(X.rows(random_samples), 0);
    sample = sum.rows(pathway_sizes - 1);
    results.elem(find(sample > reference_scores)) += 1;
  }
  arma::mat pvals(pathway_sizes.size(), cell_count);
  pvals = (results + 1) / (sampling_size + 1);
  return pvals;
}
