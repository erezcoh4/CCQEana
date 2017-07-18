// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME CCQEana_CCQEPackageDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "AnalyseEvents.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *AnalyseEvents_Dictionary();
   static void AnalyseEvents_TClassManip(TClass*);
   static void *new_AnalyseEvents(void *p = 0);
   static void *newArray_AnalyseEvents(Long_t size, void *p);
   static void delete_AnalyseEvents(void *p);
   static void deleteArray_AnalyseEvents(void *p);
   static void destruct_AnalyseEvents(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AnalyseEvents*)
   {
      ::AnalyseEvents *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AnalyseEvents));
      static ::ROOT::TGenericClassInfo 
         instance("AnalyseEvents", "AnalyseEvents.h", 27,
                  typeid(::AnalyseEvents), DefineBehavior(ptr, ptr),
                  &AnalyseEvents_Dictionary, isa_proxy, 4,
                  sizeof(::AnalyseEvents) );
      instance.SetNew(&new_AnalyseEvents);
      instance.SetNewArray(&newArray_AnalyseEvents);
      instance.SetDelete(&delete_AnalyseEvents);
      instance.SetDeleteArray(&deleteArray_AnalyseEvents);
      instance.SetDestructor(&destruct_AnalyseEvents);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AnalyseEvents*)
   {
      return GenerateInitInstanceLocal((::AnalyseEvents*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::AnalyseEvents*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AnalyseEvents_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AnalyseEvents*)0x0)->GetClass();
      AnalyseEvents_TClassManip(theClass);
   return theClass;
   }

   static void AnalyseEvents_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_AnalyseEvents(void *p) {
      return  p ? new(p) ::AnalyseEvents : new ::AnalyseEvents;
   }
   static void *newArray_AnalyseEvents(Long_t nElements, void *p) {
      return p ? new(p) ::AnalyseEvents[nElements] : new ::AnalyseEvents[nElements];
   }
   // Wrapper around operator delete
   static void delete_AnalyseEvents(void *p) {
      delete ((::AnalyseEvents*)p);
   }
   static void deleteArray_AnalyseEvents(void *p) {
      delete [] ((::AnalyseEvents*)p);
   }
   static void destruct_AnalyseEvents(void *p) {
      typedef ::AnalyseEvents current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AnalyseEvents

namespace ROOT {
   static TClass *vectorlEhitgR_Dictionary();
   static void vectorlEhitgR_TClassManip(TClass*);
   static void *new_vectorlEhitgR(void *p = 0);
   static void *newArray_vectorlEhitgR(Long_t size, void *p);
   static void delete_vectorlEhitgR(void *p);
   static void deleteArray_vectorlEhitgR(void *p);
   static void destruct_vectorlEhitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<hit>*)
   {
      vector<hit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<hit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<hit>", -2, "vector", 457,
                  typeid(vector<hit>), DefineBehavior(ptr, ptr),
                  &vectorlEhitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<hit>) );
      instance.SetNew(&new_vectorlEhitgR);
      instance.SetNewArray(&newArray_vectorlEhitgR);
      instance.SetDelete(&delete_vectorlEhitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEhitgR);
      instance.SetDestructor(&destruct_vectorlEhitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<hit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<hit>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEhitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<hit>*)0x0)->GetClass();
      vectorlEhitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEhitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEhitgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<hit> : new vector<hit>;
   }
   static void *newArray_vectorlEhitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<hit>[nElements] : new vector<hit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEhitgR(void *p) {
      delete ((vector<hit>*)p);
   }
   static void deleteArray_vectorlEhitgR(void *p) {
      delete [] ((vector<hit>*)p);
   }
   static void destruct_vectorlEhitgR(void *p) {
      typedef vector<hit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<hit>

namespace ROOT {
   static TClass *vectorlEPandoraNuTrackgR_Dictionary();
   static void vectorlEPandoraNuTrackgR_TClassManip(TClass*);
   static void *new_vectorlEPandoraNuTrackgR(void *p = 0);
   static void *newArray_vectorlEPandoraNuTrackgR(Long_t size, void *p);
   static void delete_vectorlEPandoraNuTrackgR(void *p);
   static void deleteArray_vectorlEPandoraNuTrackgR(void *p);
   static void destruct_vectorlEPandoraNuTrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<PandoraNuTrack>*)
   {
      vector<PandoraNuTrack> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<PandoraNuTrack>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<PandoraNuTrack>", -2, "vector", 457,
                  typeid(vector<PandoraNuTrack>), DefineBehavior(ptr, ptr),
                  &vectorlEPandoraNuTrackgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<PandoraNuTrack>) );
      instance.SetNew(&new_vectorlEPandoraNuTrackgR);
      instance.SetNewArray(&newArray_vectorlEPandoraNuTrackgR);
      instance.SetDelete(&delete_vectorlEPandoraNuTrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPandoraNuTrackgR);
      instance.SetDestructor(&destruct_vectorlEPandoraNuTrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<PandoraNuTrack> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<PandoraNuTrack>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPandoraNuTrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<PandoraNuTrack>*)0x0)->GetClass();
      vectorlEPandoraNuTrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPandoraNuTrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPandoraNuTrackgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<PandoraNuTrack> : new vector<PandoraNuTrack>;
   }
   static void *newArray_vectorlEPandoraNuTrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<PandoraNuTrack>[nElements] : new vector<PandoraNuTrack>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPandoraNuTrackgR(void *p) {
      delete ((vector<PandoraNuTrack>*)p);
   }
   static void deleteArray_vectorlEPandoraNuTrackgR(void *p) {
      delete [] ((vector<PandoraNuTrack>*)p);
   }
   static void destruct_vectorlEPandoraNuTrackgR(void *p) {
      typedef vector<PandoraNuTrack> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<PandoraNuTrack>

namespace {
  void TriggerDictionaryInitialization_libCCQEana_CCQEPackage_Impl() {
    static const char* headers[] = {
"AnalyseEvents.h",
0
    };
    static const char* includePaths[] = {
"/Users/erezcohen/larlite/UserDev/mySoftware",
"/Users/erezcohen/larlite/UserDev/AnalysisTreesInformation/AnaTreesPackage",
"/Users/erezcohen/larlite/UserDev/MyLarLite/MyPackage",
"/Users/erezcohen/larlite/UserDev/LarLite/MyPackage",
"/Users/erezcohen/larlite/UserDev/BasicTool/GeoAlgo",
"/Users/erezcohen/larlite/core",
"/Users/erezcohen/root6/root-6.04.10/include",
"/Users/erezcohen/larlite/UserDev/CCQEana/CCQEPackage/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$AnalyseEvents.h")))  AnalyseEvents;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "AnalyseEvents.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"AnalyseEvents", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libCCQEana_CCQEPackage",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libCCQEana_CCQEPackage_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libCCQEana_CCQEPackage_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libCCQEana_CCQEPackage() {
  TriggerDictionaryInitialization_libCCQEana_CCQEPackage_Impl();
}
