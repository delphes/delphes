#ifndef VertexFinderDA4D_h
#define VertexFinderDA4D_h

/** \class VertexFinderDA4D
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  Author O. Cerri
 *
 */

#include "classes/DelphesModule.h"

#include <vector>
#include <iostream>

class TObjArray;
class TIterator;
class Candidate;

class VertexFinderDA4D: public DelphesModule
{
  public:

    class tracks_t
    {
      public:
        std::vector<double> z;      // z-coordinate at point of closest approach to the beamline [mm]
        std::vector<double> r;      // z-coordinate at point of closest approach to the beamline [mm]

        std::vector<double> t;      // t-coordinate at point of closest approach to the beamline  [ps]

        std::vector<double> dz2_o;    //1 over the squared error of z(pca)
        std::vector<double> dz_o;
        std::vector<double> dt2_o;    //1 over the squared error of t(pca)
        std::vector<double> dt_o;

        std::vector<Candidate*> tt; // a pointer to the Candidate Track
        std::vector<double> Z;      // Normalization of P(y|x) = p_k * exp(-beta*E_ik)/Z_i

        std::vector<double> w;     // track weight

        // std::vector<double> pt;
        // std::vector<double> pz;

        std::vector<int> PID;

        double sum_w_o_dz2 = 0;
        double sum_w_o_dt2 = 0;
        double sum_w = 0;

        tracks_t(){}
        ~tracks_t(){}

        void addItem( double new_z, double new_t, double new_dz2_o, double new_dt2_o, Candidate * new_tt, double new_w, int new_PID)
        {
          z.push_back( new_z );
          t.push_back( new_t );
          dz2_o.push_back( new_dz2_o );
          dz_o.push_back( sqrt(new_dz2_o) );
          dt2_o.push_back( new_dt2_o );
          dt_o.push_back( sqrt(new_dt2_o) );
          tt.push_back( new_tt );

          w.push_back( new_w );
          Z.push_back( 1.0 );

          PID.push_back(new_PID);
        }

        unsigned int getSize()
        {
          return z.size();
        }

    };

    class vertex_t
    {
      public:
        std::vector<double> z;  //           z coordinate
        std::vector<double> t;  //           t coordinate
        std::vector<double> pk; //          vertex weight for "constrained" clustering

        // Elements of the covariance matrix of the posterior distribution
        std::vector<double> szz, stt, stz;
        std::vector<double> beta_c;

        double ZSize;
        double TSize;

        // Used in the tracks-cluster assignment
        double sum_pt;
        std::vector<unsigned int> i_tracks;


        vertex_t(){}
        ~vertex_t(){}

        void addItem( double new_z, double new_t, double new_pk)
        {
          z.push_back( new_z);
          t.push_back( new_t);
          pk.push_back( new_pk);

          szz.push_back(0);
          stt.push_back(0);
          stz.push_back(0);
          beta_c.push_back(0);
        }

        void removeItem( unsigned int i )
        {
          z.erase( z.begin() + i );
          t.erase( t.begin() + i );
          pk.erase( pk.begin() + i );

          szz.erase(szz.begin() + i);
          stt.erase(stt.begin() + i);
          stz.erase(stz.begin() + i);
          beta_c.erase(beta_c.begin() + i);
        }

        void mergeItems(unsigned int k1, unsigned int k2)
        {
          double new_pk = pk[k1] + pk[k2];
          double new_z = (z[k1]*pk[k1] + z[k2]*pk[k2])/new_pk;
          double new_t = (t[k1]*pk[k1] + t[k2]*pk[k2])/new_pk;

          removeItem(k2);
          z[k1] = new_z;
          t[k1] = new_t;
          pk[k1] = new_pk;
        }

        unsigned int getSize() const
        {
          return z.size();
        }

        double DistanceSquare(unsigned int k1, unsigned int k2)
        {
          double dz = (z[k1] - z[k2])/ZSize;
          double dt = (t[k1] - t[k2])/TSize;

          return dz*dz + dt*dt;
        }

        unsigned int NearestCluster(double t_, double z_)
        {
          unsigned int k_min = 0;
          double d2_min = 0;

          for(unsigned int k = 0; k < getSize(); k++)
          {
            double dt = (t[k] - t_)/TSize;
            double dz = (z[k] - z_)/ZSize;
            double d2 = dt*dt + dz*dz;
            if(k == 0 || d2 < d2_min)
            {
              k_min = k;
              d2_min = d2;
            }
          }

          return k_min;
        }

        std::pair<double, unsigned int> ComputeAllBeta_c(unsigned int verbose = 0)
        {
          unsigned int nv = getSize();
          unsigned int k_min = 0;
          double beta_c_min = 0;
          for(unsigned int k = 0; k < nv; k++)
          {
            if(szz[k] == 0 || stt[k] == 0)
            {
              beta_c[k] = 1000;
              if(verbose>5)
              {
                std::cout << std::endl << Form("Varianca too small! beta_c set to 1000\nVertex %d:  t=%1.2e, z=%1.2e, pk=%1.2e, szz=%1.2e, stt=%1.2e, stz=%1.2e", k, t[k], z[k], pk[k], szz[k], stt[k], stz[k]);
                // throw std::invalid_argument( "Attempting to compute beta c for uninitialized vertex" );
              }
            }

            double aux = (stt[k] - szz[k])*(stt[k] - szz[k]) + 4*stz[k]*stz[k];
            aux = 1. / (stt[k] + szz[k] + sqrt(aux));

            if(k == 0 || aux < beta_c_min)
            {
              beta_c_min = aux;
              k_min = k;
            }

            beta_c[k] = aux;
          }

          std::pair<double, unsigned int> out(beta_c_min, k_min);
          return out;
        }
    };

    VertexFinderDA4D();
    ~VertexFinderDA4D();

    void Init();
    void Process();
    void Finish();

    void clusterize(TObjArray &clusters);

    // Define the distance metric between tracks and vertices
    double Energy(double t_z, double v_z, double dz2_o, double t_t, double v_t, double dt2_o);

    // Fill the tracks structure from the input array
    void fill(tracks_t &tks);

    // Compute higher phase transition temperature
    double beta0(tracks_t & tks, vertex_t &vtx);

    // Compute the new vertexes position and mass
    double update(double beta, tracks_t &tks, vertex_t &vtx, double rho0);

    // If a vertex has beta_c lower than beta, split it
    bool split(double &beta,  vertex_t & vtx, tracks_t & tks);

    // Merge vertexes closer than declared dimensions
    bool merge(vertex_t & vtx, double d2_merge);

    // Eliminate clusters with only one significant/unique track
    bool purge(vertex_t & vtx, tracks_t & tks, double & rho0, const double beta, double min_prob, double min_trk);

    // Compute all the energies and set the partition function normalization for each track
    std::vector<double> Compute_pk_exp_mBetaE(double beta, vertex_t &vtx, tracks_t &tks, double Z_init);

    // Plot status of tracks and Vertices
    void plot_status(double beta, vertex_t &vtx, tracks_t &tks, int n_it = 0, const char* flag ="");

    // Plot status at the end of the fitting
    void plot_status_end(vertex_t &vtx, tracks_t &tks);

    // Plot Crush
    void plot_split_crush(double zn, double tn, vertex_t &vtx, tracks_t &tks, int i_vtx);


  private:

    unsigned int fVerbose;

    UInt_t fMaxIterations;
    UInt_t fMaxVertexNumber;

    Float_t fBetaMax;
    Float_t fBetaStop;
    Float_t fBetaPurge;


    Float_t fVertexZSize;
    Float_t fVertexTSize;

    Double_t fCoolingFactor;

    Double_t fDzCutOff;
    Double_t fD0CutOff;
    Double_t fDtCutOff;
    Double_t fPtMin;
    Double_t fPtMax;

    Double_t fD2UpdateLim;
    Double_t fD2Merge;
    Double_t fMuOutlayer;
    Double_t fMinTrackProb;
    Int_t fMinNTrack;

    TObjArray *fInputArray;
    TIterator *fItInputArray;

    TObjArray *fTrackOutputArray;
    TObjArray *fVertexOutputArray;

    ClassDef(VertexFinderDA4D, 1)

    //Used to keep track of number of verices in verbose mode
    std::vector<double> fEnergy_rec;
    std::vector<double> fBeta_rec;
    std::vector<double> fNvtx_rec;

    // Debug purpose only
    TObjArray *fInputGenVtx;
    TIterator *fItInputGenVtx;

    TString fFigFolderPath = "";
};

#endif
