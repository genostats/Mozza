//[[Rcpp::depends(gaston)]]
#include <Rcpp.h>
#include "sampleInt.h"
#include "Mozza.h"
#include "gaston/matrix4.h"
#include "phenotyper.h"

using namespace Rcpp;


// calibration de la composante génétique

void calibrate(int nbSimHaps, std::vector< std::vector<double> > & PrHaps, NumericVector probaDemes, NumericVector liabilityOffsetDemes,
               double length_tiles, double h2, mozza::phenotyper<IntegerVector, NumericVector> & PT) {

  int nbDemes = PrHaps.size(); // ou probaDemes.size()
  if(nbDemes == 0) stop("No demes");

  double sumG = 0, sumGG = 0;
  for(int i = 0; i < nbSimHaps; i++) {
    int d = sampleInt(probaDemes);
    mozza::zygote z = mozza::zygote(mozza::Autosomes(), PrHaps[d], length_tiles);
    double G = PT.getZygoteScore(z);
    sumG += G;
    sumGG += G*G;
  }
  // E(G) et var(G)
  double mu = sumG / (double) nbSimHaps;
  double s2 = sumGG / (double) nbSimHaps - mu*mu;

  // variance de l'offset environnemental
  double muOff = 0, mu2Off = 0;
  for(int i = 0; i < nbDemes; i++) {
    double x = probaDemes[i] * liabilityOffsetDemes[i];
    muOff += x;
    mu2Off += x * liabilityOffsetDemes[i];
  }
  double varOffset = mu2Off - muOff*muOff;
  if(varOffset > 1 - h2){
    Rcerr << "Variance of environmental offset = " << varOffset << "\n";
    stop("The variance of the environmental offset exceeds 1 - h2");
  }
Rcout << "mu Offset = " << muOff << "\n";
Rcout << "var Offset = " << varOffset << "\n";

  // s = facteur multiplicatif pour G ; sdE = variance (résiduelle) de l'environnement
  double s = sqrt(h2/s2);
  double sdE = sqrt(1 - h2 - varOffset);
  PT.setCalibration(mu + muOff/s, s, sdE); // mu + muOff/s : pour que E(G.adj) = -muOff = -E(E) 
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
              std::vector<int> & Deme, 
	      int n_target, std::vector< std::vector<double> > & PrHaps, NumericVector probaDemes, NumericVector liabilityOffsetDemes,
                      double length_tiles, double minimum_liability, double maximum_liability, 
                      mozza::phenotyper<IntegerVector, NumericVector> & PT) {

  int nbDemes = PrHaps.size(); // ou probaDemes.size()
  if(nbDemes == 0) stop("No demes");
  // int n_haps = PrHaps[0].size();

  int k = 0;
  for(int i = 0; k < n_target; i++) {
    int d = sampleInt(probaDemes);
    mozza::zygote z = mozza::zygote(mozza::Autosomes(), PrHaps[d], length_tiles);
    std::tuple<double, double, double> GGE = PT.getLiability(z);
    if(!((i+1) % 500)) Rcpp::Rcout << "\r" << k+1 << "/" << i+1 ;
    double liab = std::get<1>(GGE) + std::get<2>(GGE) + liabilityOffsetDemes[d];
    if(minimum_liability < liab && liab < maximum_liability) {
      ZYG.push_back(z);
      Graw.push_back(std::get<0>(GGE));
      G.push_back(std::get<1>(GGE));
      E.push_back(std::get<2>(GGE) + liabilityOffsetDemes[d]);
      Deme.push_back(d);
      k++;
    }
  }
  Rcpp::Rcout << "\n";
}


//[[Rcpp::export]]
List liabilitySelectedIndsProbs(IntegerVector groupSize, NumericMatrix probaHaplos, NumericVector probaDemes, 
                NumericVector liabilityOffsetDemes, 
                NumericVector minimumLiability, NumericVector maximumLiability, 
                double length_tiles, XPtr<matrix4> Haplos, IntegerVector chr, NumericVector dist, 
                IntegerVector submap, NumericVector beta, double h2, bool kinship = false, bool fraternity = false) {

  int nbGroups = groupSize.size();
  int n_haps = Haplos->ncol; // chaque haplotype = un "individu"

  if(minimumLiability.size() != nbGroups || maximumLiability.size() != nbGroups)
    stop("groupSize, minimumLiability, maximumLiability dimensions mismatch");

  int nbDemes = probaDemes.size();
  if(probaHaplos.ncol() != nbDemes || probaHaplos.nrow() != n_haps)
    stop("probaHaplos, dimensions mismatch");

  if(probaDemes.size() != nbDemes || liabilityOffsetDemes.size() != nbDemes)
    stop("probaDemes, liabilityOffsetDemes and probaHaplos dimensions mismatch");

  // chaque colonne de probaHaplos correspond à un dème
  // et donne la proba de choisir chacun des haplotype pour un individu de ce dème.
  // On crée un vecteur de std::vector<double> pour restocker ça...
  std::vector< std::vector<double> > PrHaps(nbDemes);
  for(int j = 0; j < nbDemes; j++) {
    for(int i = 0; i < n_haps; i ++) {
      PrHaps[j].push_back(probaHaplos(i,j));
    }
  }

  mozza::mappedBed<IntegerVector, NumericVector> MB(Haplos, chr, dist);
  mozza::phenotyper<IntegerVector, NumericVector> PT(MB, submap, beta);

  // calibration
  calibrate(5000, PrHaps, probaDemes, liabilityOffsetDemes, length_tiles, h2, PT);
  
  std::vector<mozza::zygote> ZYG;
  std::vector<double> Graw;
  std::vector<double> G;
  std::vector<double> E;
  std::vector<int> Deme;

  for(int i = 0; i < nbGroups; i++) {
    liabilitySelectedInds(ZYG, Graw, G, E, Deme, groupSize[i], PrHaps, probaDemes, liabilityOffsetDemes, length_tiles, minimumLiability[i], maximumLiability[i], PT); 
  }

  List L;

  Rcpp::Rcout << "Dropping to bed matrix\n";
  L["bed"] = zygote_to_bed_matrix(ZYG, MB);
  L["G.raw"] = wrap(Graw);
  L["G.adj"] = wrap(G);
  L["E"] = wrap(E);
  L["Deme"] = wrap(Deme);
  L["n"] = ZYG.size();
  if(kinship) 
    L["kinship"] = kinship_matrix(ZYG);
  if(fraternity) 
    L["fraternity"] = fraternity_matrix(ZYG);
  return L;
}

