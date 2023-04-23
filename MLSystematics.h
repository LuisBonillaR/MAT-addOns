#ifndef MLSYSTEMATICS_H
#define MLSYSTEMATICS_H

namespace PlotUtils {

//==============================================================================
// Interaction type Universe
//==============================================================================
template <class T>
class MLIntTypeUniverse : public T {
   private:
   double m_rez;
   
   //===========================================================================
   // Constructor
   //===========================================================================
   // chw: Chain with events
   // nsigma = {1,...,Univs}
   // rez: Number of universes plus 1 -> Univs + 1
   public:
   MLIntTypeUniverse(PlotUtils::ChainWrapper* chw, double nsigma, double rez) :
       m_rez(rez), T(chw, nsigma){}
   
    //==========================================================================
    // Machine learning predictions related variables
    //==========================================================================
   virtual std::vector<double> GetUnshiftedIntTypeSoft() const {
      return T::GetIntTypeSoft();
   }
    
   virtual std::vector<double> GetIntTypeSoft() const override {
      std::vector<double> softmax = T::GetIntTypeSoft();
      int maxProbIndex;
      double maxProb;
      
      if (IsMLShiftedIntType()) {
         maxProbIndex = T::GetMLIntType();
         maxProb = softmax[maxProbIndex];
         softmax[maxProbIndex] = 0.0;
         
         for (unsigned int i = 0; i < softmax.size(); ++i)
            softmax[i] /= (1 - maxProb);
      }
      
      return softmax;
   }
   
   virtual int GetUnshiftedMLIntType() const override {
      return T::GetMLIntType();
   }
   
   virtual int GetMLIntType() const override {
      double nsigma = T::m_nsigma;
      double sigRezRatio = nsigma / m_rez;
      std::vector<double> probs = T::GetIntTypeSoft();
      double nClasses = (double)probs.size();
      int maxProbIndex = T::GetMLIntType();
       
      double maxProb = probs[maxProbIndex];
      int secondMaxProbIndex;
      double secondMaxProb;
      double threshold;
      double secondThreshold;
      
      probs[maxProbIndex] = 0.0;
      secondMaxProbIndex =
         std::max_element(probs.begin(), probs.end()) - probs.begin();
      
      secondMaxProb = probs[secondMaxProbIndex] / (1 - maxProb);
      threshold = (1 + sigRezRatio * (nClasses - 1)) / nClasses;
      secondThreshold =
         (1 + sigRezRatio * (nClasses - 2)) / (nClasses - 1);
      
      // We need to take care of maxProb < 1/2 cases
      double maxAvailProb = maxProb / (1 - maxProb);
      while (secondThreshold > maxAvailProb && nsigma > 0) {
         sigRezRatio = --nsigma / m_rez;
         secondThreshold =
            (1 + sigRezRatio * (nClasses - 2)) / (nClasses - 1);
      }
      
      if (maxProb < threshold && secondMaxProb > secondThreshold)
         return secondMaxProbIndex;
      
      return maxProbIndex;
   }
   
   virtual bool IsMLShiftedIntType() const {
      return GetMLIntType() != T::GetMLIntType();
   }
   
   //===========================================================================
   // Universe information
   //===========================================================================
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
    
    //==========================================================================
    // Constructor
    //==========================================================================
   public:
   MLMultUniverse(PlotUtils::ChainWrapper* chw, double nsigma, double rez) :
      m_rez(rez), T(chw, nsigma) {}
    
   //===========================================================================
   // Machine learning predictions related variables
   //===========================================================================
   virtual std::vector<double> GetUnshiftedMultSoft() const {
      return T::GetMultSoft();
   }
    
   virtual std::vector<double> GetMultSoft() const override {
      std::vector<double> softmax = T::GetMultSoft();
      int maxProbIndex;
      double maxProb;
      
      if (IsMLShiftedMult()) {
         maxProbIndex = (int)T::GetMLHadMult();
         maxProb = softmax[maxProbIndex];
         softmax[maxProbIndex] = 0.0;
         
         for (unsigned int i = 0; i < softmax.size(); ++i)
            softmax[i] /= (1 - maxProb);
      }
      
      return softmax;
   }
   
   virtual int GetUnshiftedMLHadMult() const override {
       return T::GetMLHadMult();
   }
   
   virtual double GetMLHadMult() const override {
      double nsigma = T::m_nsigma;
      double sigRezRatio = nsigma / m_rez;
      std::vector<double> probs = T::GetMultSoft();
      double nClasses = (double)probs.size();
      int maxProbIndex = (int)T::GetMLHadMult();
     
      double maxProb = probs[maxProbIndex];
      int secondMaxProbIndex;
      double secondMaxProb;
      double threshold;
      double secondThreshold;
       
      probs[maxProbIndex] = 0.0;
      secondMaxProbIndex =
         std::max_element(probs.begin(), probs.end()) - probs.begin();
       
      secondMaxProb = probs[secondMaxProbIndex] / (1 - maxProb);
      threshold = (1 + sigRezRatio * (nClasses - 1)) / nClasses;
      secondThreshold =
         (1 + sigRezRatio * (nClasses - 2)) / (nClasses - 1);
      
      // We need to take care of maxProb < 1/2 cases
      double maxAvailProb = maxProb / (1 - maxProb);
      while (secondThreshold > maxAvailProb && nsigma > 0) {
         sigRezRatio = --nsigma / m_rez;
         secondThreshold =
            (1 + sigRezRatio * (nClasses - 2)) / (nClasses - 1);
      }
      
      if (maxProb < threshold && secondMaxProb > secondThreshold)
         return secondMaxProbIndex;
      
      return (double)maxProbIndex;
   }
   
   virtual bool IsMLShiftedMult() const {
      return GetMLHadMult() != T::GetMLHadMult();
   }
   
   //===========================================================================
   // Universe information
   //===========================================================================
   virtual std::string ShortName() const {return "ML_Reco_FSI_Mult";}
   virtual std::string LatexName() const {return "ML Reco. FSI Mult";}
   virtual bool IsVerticalOnly() const {return false;}
};


template <class T>
std::map<std::string, std::vector<T*>>
   GetMLIntTypeSystMap(typename T::config_t chw, int univs) {
   std::map<std::string, std::vector<T*>> ret;
   int rez = univs + 1;
   for (int sigma = 1; sigma < rez; ++sigma)
      ret["ML_Int_Type"].push_back(new
         PlotUtils::MLIntTypeUniverse<T>(chw, (double)sigma, (double)rez));
   
   return ret;
}


template <class T>
std::map<std::string, std::vector<T*>>
   GetMLMultMap(typename T::config_t chw, int univs) {
   std::map<std::string, std::vector<T*>> ret;
   int rez = univs + 1;
   for (int sigma = 1; sigma < rez; ++sigma)
      ret["ML_multiplicity"].push_back(new
         PlotUtils::MLMultUniverse<T>(chw, (double)sigma, (double)rez));
   
   return ret;
}

}

#endif
