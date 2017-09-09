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
#include "GENIEinteraction.h"
#include "PandoraNuTrack.h"
#include "box.h"
#include "flash.h"
#include "hit.h"
#include "pairVertex.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *hit_Dictionary();
   static void hit_TClassManip(TClass*);
   static void *new_hit(void *p = 0);
   static void *newArray_hit(Long_t size, void *p);
   static void delete_hit(void *p);
   static void deleteArray_hit(void *p);
   static void destruct_hit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::hit*)
   {
      ::hit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::hit));
      static ::ROOT::TGenericClassInfo 
         instance("hit", "hit.h", 23,
                  typeid(::hit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &hit_Dictionary, isa_proxy, 4,
                  sizeof(::hit) );
      instance.SetNew(&new_hit);
      instance.SetNewArray(&newArray_hit);
      instance.SetDelete(&delete_hit);
      instance.SetDeleteArray(&deleteArray_hit);
      instance.SetDestructor(&destruct_hit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::hit*)
   {
      return GenerateInitInstanceLocal((::hit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::hit*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *hit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::hit*)0x0)->GetClass();
      hit_TClassManip(theClass);
   return theClass;
   }

   static void hit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *box_Dictionary();
   static void box_TClassManip(TClass*);
   static void *new_box(void *p = 0);
   static void *newArray_box(Long_t size, void *p);
   static void delete_box(void *p);
   static void deleteArray_box(void *p);
   static void destruct_box(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::box*)
   {
      ::box *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::box));
      static ::ROOT::TGenericClassInfo 
         instance("box", "box.h", 22,
                  typeid(::box), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &box_Dictionary, isa_proxy, 4,
                  sizeof(::box) );
      instance.SetNew(&new_box);
      instance.SetNewArray(&newArray_box);
      instance.SetDelete(&delete_box);
      instance.SetDeleteArray(&deleteArray_box);
      instance.SetDestructor(&destruct_box);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::box*)
   {
      return GenerateInitInstanceLocal((::box*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::box*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *box_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::box*)0x0)->GetClass();
      box_TClassManip(theClass);
   return theClass;
   }

   static void box_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *PandoraNuTrack_Dictionary();
   static void PandoraNuTrack_TClassManip(TClass*);
   static void *new_PandoraNuTrack(void *p = 0);
   static void *newArray_PandoraNuTrack(Long_t size, void *p);
   static void delete_PandoraNuTrack(void *p);
   static void deleteArray_PandoraNuTrack(void *p);
   static void destruct_PandoraNuTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PandoraNuTrack*)
   {
      ::PandoraNuTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PandoraNuTrack));
      static ::ROOT::TGenericClassInfo 
         instance("PandoraNuTrack", "PandoraNuTrack.h", 54,
                  typeid(::PandoraNuTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PandoraNuTrack_Dictionary, isa_proxy, 4,
                  sizeof(::PandoraNuTrack) );
      instance.SetNew(&new_PandoraNuTrack);
      instance.SetNewArray(&newArray_PandoraNuTrack);
      instance.SetDelete(&delete_PandoraNuTrack);
      instance.SetDeleteArray(&deleteArray_PandoraNuTrack);
      instance.SetDestructor(&destruct_PandoraNuTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PandoraNuTrack*)
   {
      return GenerateInitInstanceLocal((::PandoraNuTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PandoraNuTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PandoraNuTrack_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::PandoraNuTrack*)0x0)->GetClass();
      PandoraNuTrack_TClassManip(theClass);
   return theClass;
   }

   static void PandoraNuTrack_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *pairVertex_Dictionary();
   static void pairVertex_TClassManip(TClass*);
   static void *new_pairVertex(void *p = 0);
   static void *newArray_pairVertex(Long_t size, void *p);
   static void delete_pairVertex(void *p);
   static void deleteArray_pairVertex(void *p);
   static void destruct_pairVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::pairVertex*)
   {
      ::pairVertex *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::pairVertex));
      static ::ROOT::TGenericClassInfo 
         instance("pairVertex", "pairVertex.h", 35,
                  typeid(::pairVertex), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pairVertex_Dictionary, isa_proxy, 4,
                  sizeof(::pairVertex) );
      instance.SetNew(&new_pairVertex);
      instance.SetNewArray(&newArray_pairVertex);
      instance.SetDelete(&delete_pairVertex);
      instance.SetDeleteArray(&deleteArray_pairVertex);
      instance.SetDestructor(&destruct_pairVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::pairVertex*)
   {
      return GenerateInitInstanceLocal((::pairVertex*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::pairVertex*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pairVertex_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::pairVertex*)0x0)->GetClass();
      pairVertex_TClassManip(theClass);
   return theClass;
   }

   static void pairVertex_TClassManip(TClass* ){
   }

} // end of namespace ROOT

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
         instance("AnalyseEvents", "AnalyseEvents.h", 34,
                  typeid(::AnalyseEvents), ::ROOT::Internal::DefineBehavior(ptr, ptr),
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
   static void *new_hit(void *p) {
      return  p ? new(p) ::hit : new ::hit;
   }
   static void *newArray_hit(Long_t nElements, void *p) {
      return p ? new(p) ::hit[nElements] : new ::hit[nElements];
   }
   // Wrapper around operator delete
   static void delete_hit(void *p) {
      delete ((::hit*)p);
   }
   static void deleteArray_hit(void *p) {
      delete [] ((::hit*)p);
   }
   static void destruct_hit(void *p) {
      typedef ::hit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::hit

namespace ROOT {
   // Wrappers around operator new
   static void *new_box(void *p) {
      return  p ? new(p) ::box : new ::box;
   }
   static void *newArray_box(Long_t nElements, void *p) {
      return p ? new(p) ::box[nElements] : new ::box[nElements];
   }
   // Wrapper around operator delete
   static void delete_box(void *p) {
      delete ((::box*)p);
   }
   static void deleteArray_box(void *p) {
      delete [] ((::box*)p);
   }
   static void destruct_box(void *p) {
      typedef ::box current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::box

namespace ROOT {
   // Wrappers around operator new
   static void *new_PandoraNuTrack(void *p) {
      return  p ? new(p) ::PandoraNuTrack : new ::PandoraNuTrack;
   }
   static void *newArray_PandoraNuTrack(Long_t nElements, void *p) {
      return p ? new(p) ::PandoraNuTrack[nElements] : new ::PandoraNuTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_PandoraNuTrack(void *p) {
      delete ((::PandoraNuTrack*)p);
   }
   static void deleteArray_PandoraNuTrack(void *p) {
      delete [] ((::PandoraNuTrack*)p);
   }
   static void destruct_PandoraNuTrack(void *p) {
      typedef ::PandoraNuTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PandoraNuTrack

namespace ROOT {
   // Wrappers around operator new
   static void *new_pairVertex(void *p) {
      return  p ? new(p) ::pairVertex : new ::pairVertex;
   }
   static void *newArray_pairVertex(Long_t nElements, void *p) {
      return p ? new(p) ::pairVertex[nElements] : new ::pairVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_pairVertex(void *p) {
      delete ((::pairVertex*)p);
   }
   static void deleteArray_pairVertex(void *p) {
      delete [] ((::pairVertex*)p);
   }
   static void destruct_pairVertex(void *p) {
      typedef ::pairVertex current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::pairVertex

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
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary();
   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p);
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<float> >*)
   {
      vector<vector<float> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<float> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<float> >", -2, "vector", 477,
                  typeid(vector<vector<float> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEfloatgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<vector<float> >) );
      instance.SetNew(&new_vectorlEvectorlEfloatgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEfloatgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEfloatgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<float> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<float> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<float> >*)0x0)->GetClass();
      vectorlEvectorlEfloatgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> > : new vector<vector<float> >;
   }
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> >[nElements] : new vector<vector<float> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete ((vector<vector<float> >*)p);
   }
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete [] ((vector<vector<float> >*)p);
   }
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p) {
      typedef vector<vector<float> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<float> >

namespace ROOT {
   static TClass *vectorlEpairVertexgR_Dictionary();
   static void vectorlEpairVertexgR_TClassManip(TClass*);
   static void *new_vectorlEpairVertexgR(void *p = 0);
   static void *newArray_vectorlEpairVertexgR(Long_t size, void *p);
   static void delete_vectorlEpairVertexgR(void *p);
   static void deleteArray_vectorlEpairVertexgR(void *p);
   static void destruct_vectorlEpairVertexgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pairVertex>*)
   {
      vector<pairVertex> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pairVertex>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pairVertex>", -2, "vector", 477,
                  typeid(vector<pairVertex>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairVertexgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pairVertex>) );
      instance.SetNew(&new_vectorlEpairVertexgR);
      instance.SetNewArray(&newArray_vectorlEpairVertexgR);
      instance.SetDelete(&delete_vectorlEpairVertexgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairVertexgR);
      instance.SetDestructor(&destruct_vectorlEpairVertexgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pairVertex> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<pairVertex>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairVertexgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<pairVertex>*)0x0)->GetClass();
      vectorlEpairVertexgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairVertexgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairVertexgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pairVertex> : new vector<pairVertex>;
   }
   static void *newArray_vectorlEpairVertexgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pairVertex>[nElements] : new vector<pairVertex>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairVertexgR(void *p) {
      delete ((vector<pairVertex>*)p);
   }
   static void deleteArray_vectorlEpairVertexgR(void *p) {
      delete [] ((vector<pairVertex>*)p);
   }
   static void destruct_vectorlEpairVertexgR(void *p) {
      typedef vector<pairVertex> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<pairVertex>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 477,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

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
         instance("vector<hit>", -2, "vector", 477,
                  typeid(vector<hit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
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
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<hit> : new vector<hit>;
   }
   static void *newArray_vectorlEhitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<hit>[nElements] : new vector<hit>[nElements];
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
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 477,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

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
         instance("vector<PandoraNuTrack>", -2, "vector", 477,
                  typeid(vector<PandoraNuTrack>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
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
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<PandoraNuTrack> : new vector<PandoraNuTrack>;
   }
   static void *newArray_vectorlEPandoraNuTrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<PandoraNuTrack>[nElements] : new vector<PandoraNuTrack>[nElements];
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
"GENIEinteraction.h",
"PandoraNuTrack.h",
"box.h",
"flash.h",
"hit.h",
"pairVertex.h",
0
    };
    static const char* includePaths[] = {
"/Users/erezcohen/larlite/UserDev/mySoftware",
"/Users/erezcohen/larlite/UserDev/AnalysisTreesInformation/AnaTreesPackage",
"/Users/erezcohen/larlite/UserDev/MyLarLite/MyPackage",
"/Users/erezcohen/larlite/UserDev/LarLite/MyPackage",
"/Users/erezcohen/larlite/UserDev/BasicTool/GeoAlgo",
"/Users/erezcohen/larlite/UserDev/CCQEana/LArSoftCodes",
"/Users/erezcohen/larlite/UserDev/CCQEana/LArSoftCodes/MyObjects",
"/Users/erezcohen/larlite/core",
"/usr/local/Cellar/root6/6.06.02/include/root",
"/Users/erezcohen/larlite/UserDev/CCQEana/CCQEPackage/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libCCQEana_CCQEPackage dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$AnalyseEvents.h")))  hit;
class __attribute__((annotate("$clingAutoload$AnalyseEvents.h")))  box;
class __attribute__((annotate("$clingAutoload$AnalyseEvents.h")))  PandoraNuTrack;
class __attribute__((annotate("$clingAutoload$AnalyseEvents.h")))  pairVertex;
class __attribute__((annotate("$clingAutoload$AnalyseEvents.h")))  AnalyseEvents;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libCCQEana_CCQEPackage dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "AnalyseEvents.h"
#include "GENIEinteraction.h"
#include "PandoraNuTrack.h"
#include "box.h"
#include "flash.h"
#include "hit.h"
#include "pairVertex.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"AnalyseEvents", payloadCode, "@",
"PandoraNuTrack", payloadCode, "@",
"box", payloadCode, "@",
"hit", payloadCode, "@",
"pairVertex", payloadCode, "@",
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
