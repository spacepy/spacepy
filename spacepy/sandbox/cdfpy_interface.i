/* cdfpy_interface.i */
 %module cdfpy_interface
 %include "carrays.i"
 %include "cstring.i"

 %include "cpointer.i"
  /* %include "typemaps.i"
 %include "cmalloc.i"*/
 %include "cdata.i" 


  /* void *malloc(unsigned nbytes);
     void free(void *);   */
/* %malloc(int);
 %free(int);

 %malloc(char);
 %free(char);

 %malloc(int *, intp);
 %free(int *, intp);

 %malloc(void *, voidp);
 %free(void *, voidp);

 %malloc(char *, charp);
 %free(char *, charp);
*/
 %array_class(int, intArray);
 %array_class(char, charArray);
 %array_class(long, longArray);
 %array_class(float, floatArray);

 %array_class(short, CDF_INT2Array);
 %array_class(float, CDF_REAL4Array);
 %array_class(float, CDF_FLOATArray);
 %array_class(int, CDF_INT4Array);
 %array_class(double, CDF_EPOCHArray);
 %array_class(double, CDF_REAL8Array);
 %array_class(double, CDF_DOUBLEArray);
 %array_class(char, CDF_BYTEArray);
 %array_class(char, CDF_CHARArray);
 %array_class(unsigned char, CDF_UCHARArray);


/*
 %cdata(float, float);

 %pointer_cast(int *, unsigned int *, int_to_uint);
 %pointer_cast(void *, float *, void_to_float);

 %pointer_cast(void *, float *, void_to_CDF_REAL4);
 %pointer_cast(void *, long *, void_to_CDF_INT4);
 %pointer_cast(void *, double *, void_to_CDF_EPOCH);
 %pointer_cast(void *, double *, void_to_CDF_REAL8);
 %pointer_cast(void *, unsigned char *, void_to_CDF_BYTE);
 %pointer_cast(void *, char *, void_to_CDF_INT1);
 %pointer_cast(void *, unsigned char *, void_to_CDF_UINT1);
 %pointer_cast(void *, int *, void_to_CDF_INT2);
 %pointer_cast(void *, unsigned long *, void_to_CDF_UINT4);



 %allocators(double);
*/


 %{
 /* Put header files here or function declarations like below */
 #include "/Applications/cdf33_0-dist/include/cdf.h"
 #include <stdlib.h>

 %}
 
 %cstring_bounded_output(char *attrName, CDF_ATTR_NAME_LEN256); 
 %cstring_bounded_output(char *varName, CDF_VAR_NAME_LEN256+1);
/*  %cstring_bounded_output(char *buffer, CDF_VAR_NAME_LEN256+1); 
 %cstring_output_allocate(char *buffer, free(*$1));
*/
 %include "/Applications/cdf33_0-dist/include/cdf.h"
 %pointer_functions(int, intp);
 %pointer_functions(double, doublep);
 %pointer_functions(char, charp);
 %pointer_functions(long, longp);
 %pointer_functions(CDFid, CDFidp);
/* %newobject PyCDFgetAttr;  */
/* void PyCDFgetAttr(CDFid id, long attrNum, long entryNum, char *OUTPUT);  */
 extern CDFstatus PyCDFcreate(char *INPUT, CDFid *INOUT);
 
  %inline %{
  /* Create any sort of [size] array */
  int *int_array(int size) {
     return (int *) malloc(size*sizeof(int));
  }


 CDFstatus PyCDFcreate(char *CDFname, CDFid *id) {
    return CDFcreateCDF(CDFname, id);
   }
 CDFstatus PyCDFopen(char *CDFname, CDFid *id) {
    return CDFopen(CDFname,id);
 }
 CDFstatus PyCDFcloseCDF (CDFid id) {
    return CDFcloseCDF(id);
 }
 CDFstatus PyCDFinquire (CDFid id, long *numDims, 
     long dimSizes[CDF_MAX_DIMS], long *encoding, 
     long *majority, long *maxrRec, long *numrVars, 
     long *maxzRec, long *numzVars, long *numAttrs) {
    return CDFinquireCDF(id, numDims, dimSizes, encoding, 
        majority, maxrRec, numrVars, maxzRec, numzVars, numAttrs);
 }
 CDFstatus PyCDFinquireCDF (CDFid id, long *numDims, 
     long dimSizes[CDF_MAX_DIMS], long *encoding, 
     long *majority, long *maxrRec, long *numrVars, 
     long *maxzRec, long *numzVars, long *numAttrs) {
    return CDFinquireCDF(id, numDims, dimSizes, encoding, 
        majority, maxrRec, numrVars, maxzRec, numzVars, numAttrs);
 }
 CDFstatus PyCDFgetAttrName (CDFid id, long attrNum,
     char *attrName) {
    return CDFgetAttrName(id, attrNum, attrName);
 }
 CDFstatus PyCDFinquireAttr(CDFid id, long attrNum, char *attrName, 
			    long *attrScope, long *maxgEntry, 
			    long *maxrEntry, long *maxzEntry) {
   return CDFinquireAttr(id, attrNum, attrName, attrScope, maxgEntry, 
			 maxrEntry, maxzEntry);
 }
 CDFstatus PyCDFinquirezVar (CDFid id, long varNum, char *varName, 
     long *dataType, long *numElements, long *numDims, long dimSizes[CDF_MAX_DIMS],
         long *recVary, long dimVarys[CDF_MAX_DIMS]) {
    return CDFinquirezVar(id, varNum, varName, dataType, numElements, numDims, 
        dimSizes, recVary, dimVarys);
  }
 CDFstatus PyCDFinquirerVar (CDFid id, long varNum, char *varName,
     long *dataType, long *numElements, long *recVary, 
         long dimVarys[CDF_MAX_DIMS]){
    return CDFvarInquire(id, varNum, varName, dataType, numElements, recVary, 
        dimVarys);
  }

 CDFstatus PyCDFinquireAttrgEntry (CDFid id, long attrNum, long entryNum, long *dataType, long *numElements){
   return CDFinquireAttrgEntry(id, attrNum, entryNum, dataType, numElements);
 }

 CDFstatus PyCDFgetAttrgEntry (CDFid id, long attrNum, long entryNum, char *value) {
   /* value = (char *) malloc ((numElements+1)*sizeof(char)); */
   CDFstatus tmp = CDFgetAttrgEntry(id, attrNum, entryNum, (void*)value);
   /* printf("value=%s\n", value);
      printf("status=%ld\n", tmp); */

   return tmp;
 }



 CDFstatus PyCDFcreateAttr(CDFid id, char *attrName, long attrScope, long *attrNum) {
    return CDFattrCreate(id, attrName, attrScope, attrNum);
 }
 
 CDFstatus PyCDFinquireAttrzEntry(CDFid id, long attrNum, long entryNum, long *dataType, long *numElements) {
   return CDFinquireAttrzEntry(id, attrNum, entryNum, dataType, numElements);
 } 

 CDFstatus PyCDFgetAttrzEntryDataType(CDFid id, long attrNum, long entryNum, long *dataType) {
   return CDFgetAttrzEntryDataType(id, attrNum, entryNum, dataType);
 }

 CDFstatus PyCDFgetAttrzEntryNumElements(CDFid id, long attrNum, long entryNum, long *numElems) {
   return CDFgetAttrzEntryNumElements(id, attrNum, entryNum, numElems);
 }

 CDFstatus PyCDFgetAttrzEntry(CDFid id, long attrNum, long entryNum, void *value) {
   CDFstatus tmp = CDFgetAttrzEntry(id, attrNum, entryNum, value);
   /* printf("value=%s\n", value);
      printf("status=%ld\n", tmp); */
   return tmp;
 }

 long PyCDFgetAttrNum(CDFid id, char * attrName) {
   return CDFgetAttrNum(id, attrName);
 }

 CDFstatus PyCDFgetzVarAllocRecords(CDFid id, long varNum, long *numRecs) {
   return CDFgetzVarAllocRecords(id, varNum, numRecs);
 }

 CDFstatus PyCDFgetzVarRecordData(CDFid id, long varNum, long recNum, void * buffer) {
   float *val;
   CDFstatus tmp = CDFgetzVarRecordData(id, varNum, recNum, buffer);
   val = (float*)(buffer);
   /*  printf("buffer:%f    %f", *val, 3.14159);  */
   return tmp;

 }


 CDFstatus PyCDFhyperGetzVarData(CDFid id, long varNum, long recStart, long recCount, long recInterval, 
				 long *indices, long *counts, long *intervals, void *buffer) {
   return CDFhyperGetzVarData(id, varNum, recStart, recCount, recInterval, indices, counts, intervals, buffer);
   
 }


 CDFstatus PyCDFgetCompression(CDFid id, long *compressionType, long *compressionParams, long *compressionPercentage) {
   return CDFgetCompression(id, *compressionType, *compressionParams, *compressionPercentage);
 }



 %} 
