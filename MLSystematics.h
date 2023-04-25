#ifndef MLSYSTEMATICS_H
#define MLSYSTEMATICS_H

namespace PlotUtils {


// Helper function:
// sigma: number of current universe
// rez: total number of universes + 1
// nClasses: number of classes
// probs: vector of confidence (softmax)
// maxProbIndex: reconstructed class. Will be modified after running
// This function will be used for the calculation of alternate unverses, so
// the maxProbIndex will likely change
bool CalcMLClassShift(double sigma,
                      double rez,
                      double nClasses,
                      std::vector<double> probs,
                      int& maxProbIndex)  {
  
  bool hasShift = false;
  
  maxProbIndex =
    std::max_element(probs.begin(), probs.end()) - probs.begin();
  double maxProb = probs[maxProbIndex];
  
  int secondMaxProbIndex;
  double secondMaxProb;
  
  double threshold;
  double secondThreshold;
  double sigRezRatio = sigma / rez;
  
  threshold = (1 + sigRezRatio * (nClasses - 1)) / nClasses;
  probs[maxProbIndex] = 0.0;
  secondMaxProbIndex =
    std::max_element(probs.begin(), probs.end()) - probs.begin();
  if (nClasses == 2) {
    secondMaxProb = 2 * probs[secondMaxProbIndex];
    secondThreshold = sigRezRatio;
  }
  else {
    secondMaxProb = probs[secondMaxProbIndex] / (1 - maxProb);
    secondThreshold = (1 + sigRezRatio * (nClasses - 2)) / (nClasses - 1);
  }
  
  // We need to take care of maxProb < 1/2 cases
  double maxAvailProb = maxProb / (1 - maxProb);
  while (secondThreshold > maxAvailProb && sigma > 0) {
    sigRezRatio = --sigma / rez;
    secondThreshold =
      (1 + sigRezRatio * (nClasses - 2)) / (nClasses - 1);
  }
  
  if (maxProb < threshold && secondMaxProb > secondThreshold) {
    maxProbIndex = secondMaxProbIndex;
    hasShift = true;
  }
  
  return hasShift;
}



//==============================================================================
// Interaction type Universe
//==============================================================================
template <class T>
class MLIntTypeUniverse : public T {
  private:
  double m_rez;
  int m_depth;
  mutable int m_lastEntry = -1;
  mutable int m_intType = -1;
  mutable double m_nClasses = 3;
  mutable std::vector<double> m_conf;
  
  typedef PlotUtils::ChainWrapper ChW;
   
  //============================================================================
  // Constructor
  //============================================================================
  // Inputs:
  // chw: Chain with events
  // nsigma = {1,...,Univs} current universe
  // univs: Total number of universes
  // depth: maximum number of shifts. depth = 1 allows shift to second most
  //   likely class, depth = 2 allows shit to third most likely class, etc.
  //   depth = -1 allows all shifts.
  // Internal parameters:
  //   m_rez: "resolution" of the effective confidence interval scanning
  public:
  MLIntTypeUniverse(ChW* chw, double nsigma, double univs, int depth = -1) :
    m_rez(univs+1), m_depth(depth), T(chw, nsigma) {}
  
  //============================================================================
  // Machine learning predictions related variables
  //============================================================================
  virtual std::vector<double> GetUnshiftedIntTypeSoft() const {
    return T::GetIntTypeSoft();
  }
  
  
  virtual std::vector<double> GetIntTypeSoft() const override {
    CalcMLIntTypeVars();
    
    return m_conf;
  }
  
  
  virtual int GetUnshiftedNMLIntTypeClasses() const override{
    return T::GetNMLIntTypeClasses();
  }
  
  
  virtual int GetNMLIntTypeClasses() const override {
    return m_nClasses;
  }
  
  
  virtual int GetUnshiftedMLIntType() const override {
    return T::GetMLIntType();
  }
  
  
  virtual int GetMLIntType() const override {
    CalcMLIntTypeVars();
    
    return m_intType;
  }
  
  
  virtual bool CalcMLIntTypeVars() const {
    int currentEntry = T::GetEntry();
    // Runs only if it has not run for the current entry
    if (m_lastEntry != currentEntry) {
      m_lastEntry = currentEntry;
      
      double sigma = T::m_nsigma;
      double rez = m_rez;
      double nClasses = (double)T::GetNMLIntTypeClasses();
      std::vector<double> probs = T::GetIntTypeSoft();
      int maxProbIndex = T::GetMLIntType();
      double maxProb = probs[maxProbIndex];
      
      double iRez = rez;
      double iSigma = sigma;
      bool hasShift = true;
      
      int currentDepth = 0;
      while (hasShift && nClasses > 1 && currentDepth != m_depth) {
        ++currentDepth;
        double univs = 0;
        hasShift =
          CalcMLClassShift(iSigma, iRez, nClasses, probs, maxProbIndex);
        if (hasShift) {
          univs = 0;
          int deleteIndex =
            std::max_element(probs.begin(), probs.end()) - probs.begin();
          double deleteProb = probs[deleteIndex];
          
          // Count all universes whit shift
          if (currentDepth != m_depth) {
            int dummyIndex;
            for (double jSigma = 1; jSigma < iRez; ++jSigma) {
              bool addUniv =
                CalcMLClassShift(jSigma, iRez, nClasses, probs, dummyIndex);
              if (addUniv) {
                ++univs;
                if (jSigma == iSigma) iSigma = univs;
              }
            }
            
            iRez = univs + 1;
          } // End shifted universes count
          
          // do not change probabilities when there is only two classes!
          if (nClasses > 2) {
          probs[deleteIndex] = 0.0;
          for (int index = 0; index < probs.size(); ++index)
            probs[index] /= 1 - deleteProb;
          maxProb = probs[maxProbIndex];
          }
          --nClasses;
        }
      }
      
      m_nClasses = nClasses;
      m_intType = maxProbIndex;
      m_conf = probs;
    }
    
    return true;
  }
  
  
  virtual bool IsMLShiftedIntType() const {
    return GetMLIntType() != T::GetMLIntType();
  }
  
  
  //============================================================================
  // Universe information
  //============================================================================
  virtual std::string ShortName() const { return "ML_Int_Type"; }
  virtual std::string LatexName() const { return "ML Reco. Int. Type"; }
  virtual bool IsVerticalOnly() const { return false; }
};



//==============================================================================
// Multiplicity Universe
//==============================================================================
template <class T>
class MLMultUniverse : public T {
  private:
  double m_rez;
  int m_depth;
  mutable int m_lastEntry = -1;
  mutable int m_hadMult = -1;
  mutable double m_nClasses = 3;
  mutable std::vector<double> m_conf;
  
  typedef PlotUtils::ChainWrapper ChW;
   
  //============================================================================
  // Constructor
  //============================================================================
  // Inputs:
  // chw: Chain with events
  // nsigma = {1,...,Univs} current universe
  // univs: Total number of universes
  // depth: maximum number of shifts. depth = 1 allows shift to second most
  //   likely class, depth = 2 allows shit to third most likely class, etc.
  //   depth = -1 allows all shifts.
  // Internal parameters:
  //   m_rez: "resolution" of the effective confidence interval scanning
  //============================================================================
  // Constructor
  //============================================================================
  public:
  MLMultUniverse(ChW* chw, double nsigma, double univs, int depth = -1) :
    m_rez(univs+1), m_depth(depth), T(chw, nsigma) {}
   
  //============================================================================
  // Machine learning predictions related variables
  //============================================================================
  virtual std::vector<double> GetUnshiftedHadMultSoft() const {
    return T::GetHadMultSoft();
  }
  
  
  virtual std::vector<double> GetHadMultSoft() const override {
    CalcMLHadMultVars();
    
    return m_conf;
  }
  
  
  virtual int GetUnshiftedNMLHadMultClasses() const override{
    return T::GetNMLHadMultClasses();
  }
  
  
  virtual int GetNMLHadMultClasses() const override {
    return m_nClasses;
  }
  
  
  virtual int GetUnshiftedMLHadMult() const override {
    return T::GetMLHadMult();
  }
  
  
  virtual int GetMLHadMult() const override {
    CalcMLHadMultVars();
    
    return m_hadMult;
  }
  
  
  virtual bool CalcMLHadMultVars() const {
    int currentEntry = T::GetEntry();
    // Runs only if it has not run for the current entry
    if (m_lastEntry != currentEntry) {
      m_lastEntry = currentEntry;
      
      double sigma = T::m_nsigma;
      double rez = m_rez;
      double nClasses = (double)T::GetNMLHadMultClasses();
      std::vector<double> probs = T::GetHadMultSoft();
      int maxProbIndex = T::GetMLHadMult();
      double maxProb = probs[maxProbIndex];
      
      double iRez = rez;
      double iSigma = sigma;
      bool hasShift = true;
      
      int currentDepth = 0;
      while (hasShift && nClasses > 1 && currentDepth != m_depth) {
        ++currentDepth;
        double univs = 0;
        hasShift =
          CalcMLClassShift(iSigma, iRez, nClasses, probs, maxProbIndex);
        if (hasShift) {
          univs = 0;
          int deleteIndex =
            std::max_element(probs.begin(), probs.end()) - probs.begin();
          double deleteProb = probs[deleteIndex];
          
          if (currentDepth != m_depth) {
            int dummyIndex;
            for (double jSigma = 1; jSigma < iRez; ++jSigma) {
              bool addUniv =
                CalcMLClassShift(jSigma, iRez, nClasses, probs, dummyIndex);
              if (addUniv) {
                ++univs;
                if (jSigma == iSigma) iSigma = univs;
              }
            }
            
            iRez = univs + 1;
          }
          
          probs[deleteIndex] = 0.0;
          for (int index = 0; index < probs.size(); ++index)
            probs[index] /= 1 - deleteProb;
          maxProb = probs[maxProbIndex];
          --nClasses;
          
          if (currentDepth > 1 && false) {
            std::cout << "Current depth " << currentDepth << std::endl;
            std::cout << "Max depth " << m_depth << std::endl;
            std::cout << "Univs " << univs << std::endl;
            std::cout << "nClasses " << nClasses << std::endl;
            std::cout << "Prev mult " << deleteIndex << ". Next mult "
              << maxProbIndex << std::endl;
            std::cout << "Prev prob " << deleteProb << ". Next prob "
              << maxProb << "\n" << std::endl;
          }
        }
      }
      
      m_nClasses = nClasses;
      m_hadMult = maxProbIndex;
      m_conf = probs;
    }
    
    return true;
  }
  
  
  virtual bool IsMLShiftedHadMult() const {
    return GetMLHadMult() != T::GetMLHadMult();
  }
  
  //============================================================================
  // Universe information
  //============================================================================
  virtual std::string ShortName() const {return "ML_Reco_FSI_Mult";}
  virtual std::string LatexName() const {return "ML Reco. FSI Mult";}
  virtual bool IsVerticalOnly() const {return false;}
};



template <class T>
std::map<std::string, std::vector<T*>>
  GetMLIntTypeSystMap(typename T::config_t chw, int univs) {
  std::map<std::string, std::vector<T*>> ret;
  int depth = -1;
  for (double sigma = 1; sigma < (double)univs + 1; ++sigma)
    ret["ML_Int_Type"].push_back(new
      PlotUtils::MLIntTypeUniverse<T>(chw, sigma, (double)univs, depth));
  
  return ret;
}



template <class T>
std::map<std::string, std::vector<T*>>
   GetMLMultMap(typename T::config_t chw, int univs) {
   std::map<std::string, std::vector<T*>> ret;
  int depth = -1;
  for (double sigma = 1; sigma < (double)univs + 1; ++sigma)
      ret["ML_multiplicity"].push_back(new
         PlotUtils::MLMultUniverse<T>(chw, sigma, (double)univs, depth));
   
   return ret;
}

}

#endif
