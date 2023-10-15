// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LinkDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "DataAnalysis.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *DATAEVENT_Dictionary();
   static void DATAEVENT_TClassManip(TClass*);
   static void *new_DATAEVENT(void *p = 0);
   static void *newArray_DATAEVENT(Long_t size, void *p);
   static void delete_DATAEVENT(void *p);
   static void deleteArray_DATAEVENT(void *p);
   static void destruct_DATAEVENT(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DATAEVENT*)
   {
      ::DATAEVENT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DATAEVENT));
      static ::ROOT::TGenericClassInfo 
         instance("DATAEVENT", "DataAnalysis.hh", 29,
                  typeid(::DATAEVENT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DATAEVENT_Dictionary, isa_proxy, 4,
                  sizeof(::DATAEVENT) );
      instance.SetNew(&new_DATAEVENT);
      instance.SetNewArray(&newArray_DATAEVENT);
      instance.SetDelete(&delete_DATAEVENT);
      instance.SetDeleteArray(&deleteArray_DATAEVENT);
      instance.SetDestructor(&destruct_DATAEVENT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DATAEVENT*)
   {
      return GenerateInitInstanceLocal((::DATAEVENT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DATAEVENT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DATAEVENT_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DATAEVENT*)0x0)->GetClass();
      DATAEVENT_TClassManip(theClass);
   return theClass;
   }

   static void DATAEVENT_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *dets_Dictionary();
   static void dets_TClassManip(TClass*);
   static void *new_dets(void *p = 0);
   static void *newArray_dets(Long_t size, void *p);
   static void delete_dets(void *p);
   static void deleteArray_dets(void *p);
   static void destruct_dets(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::dets*)
   {
      ::dets *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::dets));
      static ::ROOT::TGenericClassInfo 
         instance("dets", "DataAnalysis.hh", 77,
                  typeid(::dets), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &dets_Dictionary, isa_proxy, 4,
                  sizeof(::dets) );
      instance.SetNew(&new_dets);
      instance.SetNewArray(&newArray_dets);
      instance.SetDelete(&delete_dets);
      instance.SetDeleteArray(&deleteArray_dets);
      instance.SetDestructor(&destruct_dets);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::dets*)
   {
      return GenerateInitInstanceLocal((::dets*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::dets*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *dets_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::dets*)0x0)->GetClass();
      dets_TClassManip(theClass);
   return theClass;
   }

   static void dets_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_DATAEVENT(void *p) {
      return  p ? new(p) ::DATAEVENT : new ::DATAEVENT;
   }
   static void *newArray_DATAEVENT(Long_t nElements, void *p) {
      return p ? new(p) ::DATAEVENT[nElements] : new ::DATAEVENT[nElements];
   }
   // Wrapper around operator delete
   static void delete_DATAEVENT(void *p) {
      delete ((::DATAEVENT*)p);
   }
   static void deleteArray_DATAEVENT(void *p) {
      delete [] ((::DATAEVENT*)p);
   }
   static void destruct_DATAEVENT(void *p) {
      typedef ::DATAEVENT current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DATAEVENT

namespace ROOT {
   // Wrappers around operator new
   static void *new_dets(void *p) {
      return  p ? new(p) ::dets : new ::dets;
   }
   static void *newArray_dets(Long_t nElements, void *p) {
      return p ? new(p) ::dets[nElements] : new ::dets[nElements];
   }
   // Wrapper around operator delete
   static void delete_dets(void *p) {
      delete ((::dets*)p);
   }
   static void deleteArray_dets(void *p) {
      delete [] ((::dets*)p);
   }
   static void destruct_dets(void *p) {
      typedef ::dets current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::dets

namespace ROOT {
   static TClass *vectorlEunsignedsPshortgR_Dictionary();
   static void vectorlEunsignedsPshortgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPshortgR(void *p = 0);
   static void *newArray_vectorlEunsignedsPshortgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPshortgR(void *p);
   static void deleteArray_vectorlEunsignedsPshortgR(void *p);
   static void destruct_vectorlEunsignedsPshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned short>*)
   {
      vector<unsigned short> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned short>", -2, "vector", 214,
                  typeid(vector<unsigned short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPshortgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned short>) );
      instance.SetNew(&new_vectorlEunsignedsPshortgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPshortgR);
      instance.SetDelete(&delete_vectorlEunsignedsPshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPshortgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned short> >()));

      ::ROOT::AddClassAlternate("vector<unsigned short>","std::vector<unsigned short, std::allocator<unsigned short> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned short>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned short>*)0x0)->GetClass();
      vectorlEunsignedsPshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPshortgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short> : new vector<unsigned short>;
   }
   static void *newArray_vectorlEunsignedsPshortgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned short>[nElements] : new vector<unsigned short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPshortgR(void *p) {
      delete ((vector<unsigned short>*)p);
   }
   static void deleteArray_vectorlEunsignedsPshortgR(void *p) {
      delete [] ((vector<unsigned short>*)p);
   }
   static void destruct_vectorlEunsignedsPshortgR(void *p) {
      typedef vector<unsigned short> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned short>

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
         instance("vector<float>", -2, "vector", 214,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

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
   static TClass *vectorlEdetsgR_Dictionary();
   static void vectorlEdetsgR_TClassManip(TClass*);
   static void *new_vectorlEdetsgR(void *p = 0);
   static void *newArray_vectorlEdetsgR(Long_t size, void *p);
   static void delete_vectorlEdetsgR(void *p);
   static void deleteArray_vectorlEdetsgR(void *p);
   static void destruct_vectorlEdetsgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<dets>*)
   {
      vector<dets> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<dets>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<dets>", -2, "vector", 214,
                  typeid(vector<dets>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdetsgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<dets>) );
      instance.SetNew(&new_vectorlEdetsgR);
      instance.SetNewArray(&newArray_vectorlEdetsgR);
      instance.SetDelete(&delete_vectorlEdetsgR);
      instance.SetDeleteArray(&deleteArray_vectorlEdetsgR);
      instance.SetDestructor(&destruct_vectorlEdetsgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<dets> >()));

      ::ROOT::AddClassAlternate("vector<dets>","std::vector<dets, std::allocator<dets> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<dets>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdetsgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<dets>*)0x0)->GetClass();
      vectorlEdetsgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdetsgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdetsgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<dets> : new vector<dets>;
   }
   static void *newArray_vectorlEdetsgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<dets>[nElements] : new vector<dets>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdetsgR(void *p) {
      delete ((vector<dets>*)p);
   }
   static void deleteArray_vectorlEdetsgR(void *p) {
      delete [] ((vector<dets>*)p);
   }
   static void destruct_vectorlEdetsgR(void *p) {
      typedef vector<dets> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<dets>

namespace ROOT {
   static TClass *vectorlEDATAEVENTgR_Dictionary();
   static void vectorlEDATAEVENTgR_TClassManip(TClass*);
   static void *new_vectorlEDATAEVENTgR(void *p = 0);
   static void *newArray_vectorlEDATAEVENTgR(Long_t size, void *p);
   static void delete_vectorlEDATAEVENTgR(void *p);
   static void deleteArray_vectorlEDATAEVENTgR(void *p);
   static void destruct_vectorlEDATAEVENTgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<DATAEVENT>*)
   {
      vector<DATAEVENT> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<DATAEVENT>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<DATAEVENT>", -2, "vector", 214,
                  typeid(vector<DATAEVENT>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDATAEVENTgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<DATAEVENT>) );
      instance.SetNew(&new_vectorlEDATAEVENTgR);
      instance.SetNewArray(&newArray_vectorlEDATAEVENTgR);
      instance.SetDelete(&delete_vectorlEDATAEVENTgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDATAEVENTgR);
      instance.SetDestructor(&destruct_vectorlEDATAEVENTgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<DATAEVENT> >()));

      ::ROOT::AddClassAlternate("vector<DATAEVENT>","std::vector<DATAEVENT, std::allocator<DATAEVENT> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<DATAEVENT>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDATAEVENTgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<DATAEVENT>*)0x0)->GetClass();
      vectorlEDATAEVENTgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDATAEVENTgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDATAEVENTgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DATAEVENT> : new vector<DATAEVENT>;
   }
   static void *newArray_vectorlEDATAEVENTgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DATAEVENT>[nElements] : new vector<DATAEVENT>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDATAEVENTgR(void *p) {
      delete ((vector<DATAEVENT>*)p);
   }
   static void deleteArray_vectorlEDATAEVENTgR(void *p) {
      delete [] ((vector<DATAEVENT>*)p);
   }
   static void destruct_vectorlEDATAEVENTgR(void *p) {
      typedef vector<DATAEVENT> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<DATAEVENT>

namespace {
  void TriggerDictionaryInitialization_LinkDict_Impl() {
    static const char* headers[] = {
"DataAnalysis.hh",
0
    };
    static const char* includePaths[] = {
"/opt/ROOT/root62006/include/",
"/wuhongyi1/zhangjizhi/LaBr3/pileupcode/0Check/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LinkDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate("$clingAutoload$DataAnalysis.hh")))  DATAEVENT;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
struct __attribute__((annotate("$clingAutoload$DataAnalysis.hh")))  dets;
typedef struct DATAEVENT DATAEVENT __attribute__((annotate("$clingAutoload$DataAnalysis.hh"))) ;
typedef struct dets dets __attribute__((annotate("$clingAutoload$DataAnalysis.hh"))) ;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LinkDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "DataAnalysis.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"DATAEVENT", payloadCode, "@",
"dets", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LinkDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LinkDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LinkDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LinkDict() {
  TriggerDictionaryInitialization_LinkDict_Impl();
}
