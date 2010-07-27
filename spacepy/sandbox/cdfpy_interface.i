/* cdfpy_interface.i */
 %module cdfpy_interface
 %include "carrays.i"
 %include "cstring.i"

 %include "cpointer.i"
  /* %include "typemaps.i" */
 %include "cmalloc.i"
 %include "cdata.i" 


 %array_functions(int, intArray);
 %array_functions(char, charArray);
 %array_functions(long, longArray);
 %array_functions(float, floatArray);

 %array_functions(short, CDF_INT2Array);
 %array_functions(float, CDF_REAL4Array);
 %array_functions(float, CDF_FLOATArray);
 %array_functions(int, CDF_INT4Array);
 %array_functions(double, CDF_EPOCHArray);
 %array_functions(double, CDF_REAL8Array);
 %array_functions(double, CDF_DOUBLEArray);
 %array_functions(char, CDF_BYTEArray);
 %array_functions(char, CDF_CHARArray);
 %array_functions(unsigned char, CDF_UCHARArray);


 %allocators(void, void);
 %allocators(void *, voidp);

 %{
 /* Put header files here or function declarations like below */
 #include "/Applications/cdf33_0-dist/include/cdf.h"
 #include "/usr/include/stdlib.h"

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



 CDFstatus PyCDFcreateAttr(CDFid id, char *attrNameIN, long attrScope, long *attrNum) {
    return CDFattrCreate(id, attrNameIN, attrScope, attrNum);
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
   return CDFgetzVarRecordData(id, varNum, recNum, buffer);
 }


 CDFstatus PyCDFhyperGetzVarData(CDFid id, long varNum, long recStart, long recCount, long recInterval, 
				 long *indices, long *counts, long *intervals, void *buffer) {
   return CDFhyperGetzVarData(id, varNum, recStart, recCount, recInterval, indices, counts, intervals, buffer);
   
 }


 CDFstatus PyCDFgetCompression(CDFid id, long *compressionType, long *compressionParams, long *compressionPercentage) {
   return CDFgetCompression(id, *compressionType, *compressionParams, *compressionPercentage);
 }

 CDFstatus PyCDFputAttrgEntry(CDFid id, long attrNum, long entryNum, long dataType, long numElements, void *value) {
   return CDFputAttrgEntry(id, attrNum, entryNum, dataType, numElements, value);
 }



 %} 
