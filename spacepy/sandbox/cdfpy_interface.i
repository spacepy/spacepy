/* cdfpy_interface.i */
 %module cdfpy_interface
 %include "cpointer.i"
 %include "typemaps.i"
 %include "cstring.i"
 
 %{
 /* Put header files here or function declarations like below */
 #include "/home/smorley/include/cdf.h"
 #include <stdlib.h>

 %}
 
 %cstring_bounded_output(char *attrName, CDF_ATTR_NAME_LEN256);
 %cstring_bounded_output(char *varName, CDF_VAR_NAME_LEN256+1);
/*  %cstring_bounded_output(char *buffer, CDF_VAR_NAME_LEN256+1); */
 %cstring_output_allocate(char *buffer, free(*$1));
 
 %include "/home/smorley/include/cdf.h"
 %pointer_functions(int, intp);
 %pointer_functions(double, doublep);
 %pointer_functions(char, charp);
 %pointer_functions(long, longp);
 %pointer_functions(CDFid, CDFidp);
 %newobject PyCDFgetAttr;
 void PyCDFgetAttr(CDFid id, long attrNum, long entryNum, char *OUTPUT);
 extern CDFstatus PyCDFcreate(char *INPUT, CDFid *INOUT);
 
  %inline %{
 CDFstatus PyCDFcreate(char *CDFname, CDFid *id) {
    return CDFcreateCDF(CDFname, id);
   }
 CDFstatus PyCDFopen(char *CDFname, CDFid *id) {
    return CDFopen(CDFname,id);
 }
 CDFstatus PyCDFclose (CDFid id) {
    return CDFclose(id);
 }
 CDFstatus PyCDFinquire (CDFid id, long *numDims, 
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
     long *attrScope, long *maxEntry) {
    return CDFattrInquire(id, attrNum, attrName, attrScope, maxEntry);
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
/* CDFstatus PyCDFgetrVar (CDFid id, ){
    return CDFvarHyperGet(id, varNum, recStart, recCount, recInterval, indices, value);
 }*/
 void PyCDFgetAttr (CDFid id, long attrNum, long entryNum, char *buffer){
    CDFstatus a, b;
    long dataType;
    long numElements;
    a = CDFattrEntryInquire(id, attrNum, entryNum, &dataType, &numElements);
    buffer = (char *) malloc (numElements+1);
    b = CDFattrGet(id, attrNum, entryNum, buffer);
    buffer[numElements] = '\0';
 }
 CDFstatus PyCDFcreateAttr (CDFid id, char *attrName, long attrScope, long *attrNum) {
    return CDFattrCreate(id, attrName, attrScope, attrNum);
 }
 
 %} 
 