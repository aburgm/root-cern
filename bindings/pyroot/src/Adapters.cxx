// Bindings
#include "PyROOT.h"
#include "Adapters.h"
#include "Utility.h"

// ROOT
#include "TInterpreter.h"
#include "TBaseClass.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TMethod.h"
#include "TFunction.h"
#include "TMethodArg.h"
#include "TList.h"
#include "TError.h"


//= TReturnTypeAdapter =======================================================
std::string PyROOT::TReturnTypeAdapter::Name( unsigned int mod ) const
{
// get the name of the return type that is being adapted
   std::string name = fName;

   if ( ! ( mod & Rflx::QUALIFIED ) )
      name = TClassEdit::CleanType( fName.c_str(), 1 );

   if ( mod & Rflx::FINAL )
      name = Utility::ResolveTypedef( name );

   return name;
}


//= TMemberAdapter ===========================================================
PyROOT::TMemberAdapter::TMemberAdapter( TMethod* meth ) : fMember( meth )
{
   /* empty */
}

//____________________________________________________________________________
PyROOT::TMemberAdapter::operator TMethod*() const
{
// cast the adapter to a TMethod* being adapted, returns 0 on failure
   return dynamic_cast< TMethod* >( const_cast< TDictionary* >( fMember ) );
}

//____________________________________________________________________________
PyROOT::TMemberAdapter::TMemberAdapter( TFunction* func ) : fMember( func )
{
   /* empty */
}

//____________________________________________________________________________
PyROOT::TMemberAdapter::operator TFunction*() const
{
// cast the adapter to a TFunction* being adapted, returns 0 on failure
   return dynamic_cast< TFunction* >( const_cast< TDictionary* >( fMember ) );
}

//____________________________________________________________________________
PyROOT::TMemberAdapter::TMemberAdapter( TDataMember* mb ) : fMember( mb )
{
   /* empty */
}

//____________________________________________________________________________
PyROOT::TMemberAdapter::operator TDataMember*() const
{
// cast the adapter to a TDataMember* being adapted, returns 0 on failure
   return dynamic_cast< TDataMember* >( const_cast< TDictionary* >( fMember ) );
}

//____________________________________________________________________________
PyROOT::TMemberAdapter::TMemberAdapter( TMethodArg* ma ) : fMember( ma )
{
   /* empty */
}

//____________________________________________________________________________
PyROOT::TMemberAdapter::operator TMethodArg*() const
{
// cast the adapter to a TMethodArg* being adapted, returns 0 on failure
   return dynamic_cast< TMethodArg* >( const_cast< TDictionary* >( fMember ) );
}

//____________________________________________________________________________
std::string PyROOT::TMemberAdapter::Name( unsigned int mod ) const
{
// Return name of the type described by fMember
   TMethodArg* arg = (TMethodArg*)*this;

   if ( arg ) {

      std::string name = arg->GetTypeNormalizedName();
      if ( mod & Rflx::QUALIFIED )
         return name;

      if ( mod & Rflx::FINAL )
         name = Utility::ResolveTypedef( name );

      return name;

   } else if ( mod & Rflx::FINAL )
      return Utility::ResolveTypedef( fMember->GetName() );

// CLING WORKAROUND -- TMethod names are fully scoped ...
   TMethod* m = (TMethod*)*this;
   if (m && m->GetClass()) {
      std::string scoped_name = m->GetName();
      std::string class_name = m->GetClass()->GetName();
      std::string::size_type pos = scoped_name.find(class_name + "::");
      if (pos == 0)      // only accept found at start
         scoped_name = scoped_name.substr(class_name.size() + 2 /* for :: */, std::string::npos);

   // CLING WORKAROUND (ROOT-5204) -- operator _Bool -> operator bool
      if ( scoped_name == "operator _Bool" )
         scoped_name = "operator bool";
   // -- CLING WORKAROUND
      return scoped_name;
   }

// CLING WORKAROUND -- should not be null, but can be due #100389
   if ( fMember != 0 )
      return fMember->GetName();
   return "<unknown>";
}

//____________________________________________________________________________
Bool_t PyROOT::TMemberAdapter::IsEnum() const
{
// test if the adapted member is of an enum type
   return fMember->Property() & kIsEnum;
}

//____________________________________________________________________________
Bool_t PyROOT::TMemberAdapter::IsPublic() const
{
// test if the adapted member represents an public (data) member
   return fMember->Property() & kIsPublic;
}

//____________________________________________________________________________
Bool_t PyROOT::TMemberAdapter::IsStatic() const
{
// test if the adapted member represents a class (data) member
   return fMember->Property() & kIsStatic;
}

//____________________________________________________________________________
size_t PyROOT::TMemberAdapter::FunctionParameterSize( Bool_t required ) const
{
// get the total number of parameters that the adapted function/method takes
   TFunction* func = (TFunction*)fMember;
   if ( ! func )
      return 0;

   if ( required == true )
      return func->GetNargs() - func->GetNargsOpt();

   return func->GetNargs();
}

//____________________________________________________________________________
PyROOT::TMemberAdapter PyROOT::TMemberAdapter::FunctionParameterAt( size_t nth ) const
{
// get the type info of the function parameter at position nth
   return (TMethodArg*)((TFunction*)fMember)->GetListOfMethodArgs()->At( nth );
}

//____________________________________________________________________________
std::string PyROOT::TMemberAdapter::FunctionParameterNameAt( size_t nth ) const
{
// get the formal name, if available, of the function parameter at position nth
   const char* name =
      ((TMethodArg*)((TFunction*)fMember)->GetListOfMethodArgs()->At( nth ))->GetName();

   if ( name )
      return name;
   return "";
}

//____________________________________________________________________________
std::string PyROOT::TMemberAdapter::FunctionParameterDefaultAt( size_t nth ) const
{
// get the default value, if available, of the function parameter at position nth
   TMethodArg* arg = (TMethodArg*)((TFunction*)fMember)->GetListOfMethodArgs()->At( nth );
   const char* def = arg->GetDefault();

   if ( ! def )
      return "";

// special case for strings: "some value" -> ""some value"
   if ( strstr( Utility::ResolveTypedef( arg->GetTypeNormalizedName() ).c_str(), "char*" ) ) {
      std::string sdef = "\"";
      sdef += def;
      sdef += "\"";
      return sdef;
   }

   return def;
}

//____________________________________________________________________________
PyROOT::TReturnTypeAdapter PyROOT::TMemberAdapter::ReturnType() const
{
// get the return type of the wrapped function/method
   return TReturnTypeAdapter( ((TFunction*)fMember)->GetReturnTypeNormalizedName() );
}

//____________________________________________________________________________
PyROOT::TScopeAdapter PyROOT::TMemberAdapter::DeclaringScope() const
{
// get the declaring scope (class) of the wrapped function/method
   TMethod* method = (TMethod*)*this;
   if ( method )
      return method->GetClass();

// happens for free-standing functions (i.e. global scope)
   return std::string( "" );
}


//= TBaseAdapter =============================================================
std::string PyROOT::TBaseAdapter::Name() const
{
// get the name of the base class that is being adapted
   return fBase->GetName();
}


//= TScopeAdapter ============================================================
PyROOT::TScopeAdapter::TScopeAdapter( TClass* klass ) : fClass( klass )
{
// wrap a class (scope)
   if ( fClass.GetClass() != 0 )
      fName = fClass->GetName();
}

//____________________________________________________________________________
PyROOT::TScopeAdapter::TScopeAdapter( const std::string& name ) :
   fClass( name.c_str() ), fName( name )
{
   /* empty */
}

PyROOT::TScopeAdapter::TScopeAdapter( const TMemberAdapter& mb ) :
      fClass( mb.Name( Rflx::SCOPED ).c_str() ),
      fName( mb.Name( Rflx::QUALIFIED | Rflx::SCOPED ) )
{
   /* empty */
}

//____________________________________________________________________________
PyROOT::TScopeAdapter PyROOT::TScopeAdapter::ByName( const std::string& name, Bool_t quiet )
{
// lookup a scope (class) by name
   Int_t oldEIL = gErrorIgnoreLevel;
   if ( quiet )
      gErrorIgnoreLevel = 3000;

   TClassRef klass( name.c_str() );
   if (klass.GetClass() && klass->GetListOfAllPublicMethods()->GetSize() == 0) {
   // sometimes I/O interferes, leading to zero methods: reload from CINT
      ClassInfo_t* cl = gInterpreter->ClassInfo_Factory( name.c_str() );
      if ( cl ) gInterpreter->SetClassInfo( klass, kTRUE );
   }

   gErrorIgnoreLevel = oldEIL;

   return klass.GetClass();
}

//____________________________________________________________________________
std::string PyROOT::TScopeAdapter::Name( unsigned int mod ) const
{
// Return name of type described by fClass
   if ( ! fClass.GetClass() || ! fClass->Property() ) {
   // fundamental types have no class, and unknown classes have no property
      std::string name = fName;

      if ( ! ( mod & Rflx::QUALIFIED ) )
         name = TClassEdit::CleanType( fName.c_str(), 1 );

      if ( mod & Rflx::FINAL )
         name = Utility::ResolveTypedef( name );

      return name;
   }

   if ( mod & Rflx::FINAL ) {
   /* The following G__ClassInfo lookup is badly broken ...
      G__ClassInfo* clInfo = (G__ClassInfo*)fClass->GetClassInfo();
      if ( mod & Rflx::SCOPED )
         return (clInfo && clInfo->Fullname()) ? clInfo->Fullname() : fClass->GetName();

   // unscoped name ...
      std::string actual = (clInfo && clInfo->Name()) ? clInfo->Name() : fClass->GetName();

   // in case of missing dictionaries, the scope won't have been stripped
      if ( ! ( clInfo && clInfo->IsValid() ) ) {
         std::string::size_type pos = actual.substr( 0, actual.find( '<' ) ).rfind( "::" );
         if ( pos != std::string::npos ) {
         // this is somewhat of a gamble, but the alternative is a guaranteed crash
            actual = actual.substr( pos + 2, std::string::npos );
         }
      }

      return actual;
   */
      return fClass->GetName();

   } else if ( ! ( mod & Rflx::SCOPED ) ) {
   /* The following G__ClassInfo lookup is badly broken ...
      G__ClassInfo* clInfo = (G__ClassInfo*)fClass->GetClassInfo();
      return (clInfo && clInfo->Name()) ? clInfo->Name() : fClass->GetName();
   */
      return fClass->GetName();
   }

   return fClass->GetName();
}

//____________________________________________________________________________
size_t PyROOT::TScopeAdapter::BaseSize() const
{
// get the total number of base classes that this class has
   if ( fClass.GetClass() && fClass->GetListOfBases() != 0 )
      return fClass->GetListOfBases()->GetSize();

   return 0;
}

//____________________________________________________________________________
PyROOT::TBaseAdapter PyROOT::TScopeAdapter::BaseAt( size_t nth ) const
{
// get the nth base of this class
   return (TBaseClass*)fClass->GetListOfBases()->At( nth );
}

//____________________________________________________________________________
size_t PyROOT::TScopeAdapter::FunctionMemberSize() const
{
// get the total number of methods that this class has
   if ( fClass.GetClass() )
      return fClass->GetListOfMethods()->GetSize();

   return 0;
}

//____________________________________________________________________________
PyROOT::TMemberAdapter PyROOT::TScopeAdapter::FunctionMemberAt( size_t nth ) const
{
// get the nth method of this class
   return (TMethod*)fClass->GetListOfMethods()->At( nth );
}

//____________________________________________________________________________
size_t PyROOT::TScopeAdapter::DataMemberSize() const
{
// get the total number of data members that this class has
   if ( fClass.GetClass() )
      return fClass->GetListOfDataMembers()->GetSize();

   return 0;
}

//____________________________________________________________________________
PyROOT::TMemberAdapter PyROOT::TScopeAdapter::DataMemberAt( size_t nth ) const
{
// get the nth data member of this class
   return (TDataMember*)fClass->GetListOfDataMembers()->At( nth );
}

//____________________________________________________________________________
PyROOT::TScopeAdapter::operator Bool_t() const
{
// check the validity of this scope (class)
   if ( fName.empty() )
      return false;

   Bool_t b = kFALSE;

   Int_t oldEIL = gErrorIgnoreLevel;
   gErrorIgnoreLevel = 3000;
   TClass* klass = TClass::GetClass( Name( Rflx::QUALIFIED | Rflx::SCOPED ).c_str() );
   if ( klass && klass->GetClassInfo() )
      b = gInterpreter->ClassInfo_IsValid( klass->GetClassInfo() );
   gErrorIgnoreLevel = oldEIL;
   return b;
}
      
//____________________________________________________________________________
Bool_t PyROOT::TScopeAdapter::IsComplete() const
{
// verify whether the dictionary of this class is fully available
   Bool_t b = kFALSE;

   TClass* klass = TClass::GetClass( Name( Rflx::SCOPED ).c_str() );
   if ( klass ) b = gInterpreter->ClassInfo_IsLoaded( klass->GetClassInfo() );
   return b;
}

//____________________________________________________________________________
Bool_t PyROOT::TScopeAdapter::IsClass() const
{
// test if this scope represents a class
   if ( fClass.GetClass() ) {
   // some inverted logic: we don't have a TClass, but a builtin will be recognized, so
   // if it is NOT a builtin, it is a class or struct (but may be missing dictionary)
      return (fClass->Property() & kIsClass) || ! (fClass->Property() & kIsFundamental);
   }

// no class can mean either is no class (i.e. builtin), or no dict but coming in
// through PyCintex/Reflex ... as a workaround, use TDataTypes that has a full
// enumeration of builtin types
   return TDataType( Name( Rflx::FINAL | Rflx::SCOPED ).c_str() ).GetType() == kOther_t;
}

//____________________________________________________________________________
Bool_t PyROOT::TScopeAdapter::IsStruct() const
{
// test if this scope represents a struct
   if ( fClass.GetClass() ) {
   // same logic as for IsClass() above ...
      return (fClass->Property() & kIsStruct) || ! (fClass->Property() & kIsFundamental);
   }

// same logic as for IsClass() above ...
   return TDataType( Name( Rflx::FINAL | Rflx::SCOPED ).c_str() ).GetType() == kOther_t;
}

//____________________________________________________________________________
Bool_t PyROOT::TScopeAdapter::IsNamespace() const
{
// test if this scope represents a namespace
   if ( fClass.GetClass() )
      return fClass->Property() & kIsNamespace;

   return kFALSE;
}

//____________________________________________________________________________
Bool_t PyROOT::TScopeAdapter::IsAbstract() const
{
// test if this scope represents an abstract class
   if ( fClass.GetClass() )
      return fClass->Property() & kIsAbstract;   // assume set only for classes

   return kFALSE;
}
