//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "mozza.h"
#include "gaston/matrix4.h"
#include "phenotyper.h"

using namespace Rcpp;


// calibration de la composante génétique

/* Une possibilité non retenue mais je garde le code sous le coude

    std::vector<double> scores = scoreAllHaplotypes<double>(MB, submap, beta);
    int n_haps = MB.nbHaps; 
    mu = std::accumulate(scores.begin(), scores.end(), 0.) / (double) n_haps;
    s2 = std::accumulate(scores.begin(), scores.end(), 0., [](double a, double b) {return a + b*b;}) / (double) n_haps - mu*mu;
*/


void calibrate(int nbSimHaps, int n_haps, double length_tiles, double h2, mozza::phenotyper<IntegerVector, NumericVector> & PT) {
  double sumG = 0, sumGG = 0;
  for(int i = 0; i < nbSimHaps; i++) {
    mozza::zygote z = mozza::zygote(mozza::human_autosomes_b37, n_haps, length_tiles);
    double G = PT.getZygoteScore(z);
    sumG += G;
    sumGG += G*G;
  }
  double mu = sumG / (double) n_haps;
  double s2 = sumGG / (double) n_haps - mu*mu;
  double s = sqrt(h2/s2);
  double sdE = sqrt(1 - h2);
  PT.setCalibration(mu, s, sdE);
}


// cette fonction fait des individus indépendants, sélectionnés pour leur liabilité (rejection sampling)
// avec haplos mosaiques avec tuiles de longueurs length_tiles parmi n_haps haplotypes
//
// Ils sont pushed back dans un vecteur de mozza::zygote
// ZYG, G, E : vector in which to push pack the zygote / les valeurs de G et E
// n_target : number of inds to generate
// n_haps : tiles will be numbered from 0 to n_haps-1
// length_tiles : in cM
// minimum_liability, maximum_liability : le critère de sélection
// PT : le phenotyper

void liabilitySelectedInds(std::vector<mozza::zygote> & ZYG, std::vector<double> & Graw, std::vector<double> & G, std::vector<double> & E, 
                      int n_target, int n_haps, double length_tiles, double minimum_liability, double maximum_liability, 
                      mozza::phenotyper<IntegerVector, NumericVector> & PT) {
  int k = 0;
  for(int i = 0; k < n_target; i++) {
    mozza::zygote z = mozza::zygote(mozza::human_autosomes_b37, n_haps, length_tiles);
    std::tuple<double, double, double> GGE = PT.getLiability(z);
    if(!((i+1) % 500)) Rcpp::Rcout << k+1 << "/" << i+1 << "\r";
    double liab = std::get<1>(GGE) + std::get<2>(GGE);
    if(minimum_liability < liab && liab < maximum_liability) {
      ZYG.push_back(z);
      Graw.push_back(std::get<0>(GGE));
      G.push_back(std::get<1>(GGE));
      E.push_back(std::get<2>(GGE));
      k++;
    }
  }
  Rcpp::Rcout << "\n";
}

//[[Rcpp::export]]
List liabilitySelectedInds(IntegerVector groupSize, NumericVector minimumLiability, NumericVector maximumLiability, 
                double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, 
                IntegerVector submap, NumericVector beta, double h2, bool kinship = false, bool fraternity = false) {

  int nbGroups = groupSize.size();
  if(minimumLiability.size() != nbGroups || maximumLiability.size() != nbGroups)
    stop("groupSize, minimumLiability, maximumLiability dimensions mismatch");

  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"
  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  mozza::phenotyper<IntegerVector, NumericVector> PT(MB, submap, beta);

  // calibration
  calibrate(1000, n_haps, length_tiles, h2, PT);
  
  std::vector<mozza::zygote> ZYG;
  std::vector<double> Graw;
  std::vector<double> G;
  std::vector<double> E;

  for(int i = 0; i < nbGroups; i++) {
    liabilitySelectedInds(ZYG, Graw, G, E, groupSize[i], n_haps, length_tiles, minimumLiability[i], maximumLiability[i], PT); 
  }

  List L;

  Rcpp::Rcout << "Dropping to bed matrix\n";
  L["bed"] = drop_to_bed_matrix(ZYG, MB);
  L["G.raw"] = wrap(Graw);
  L["G.adj"] = wrap(G);
  L["E"] = wrap(E);
  L["n"] = ZYG.size();
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  if(fraternity) 
    L["fraternity"] = fraternity_matrix(ZYG);
  return L;
}

