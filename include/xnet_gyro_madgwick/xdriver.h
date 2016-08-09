 /*
 *  file : xdriver.h
 *  author : TimXu,AllenWang
 *  note : read/write tag from controller
 */

#ifndef _H_XDRIVER
#define _H_XDRIVER

#ifdef __cplusplus
extern "C" {
#endif  /* __cplusplus */


	#define MAKEIPADDRESS(b1,b2,b3,b4) ( ((unsigned int)b1<<24)+( (unsigned int)b2<<16)+ ((unsigned int)b3<<8)+b4 )

	#define CONTROLLER_IP MAKEIPADDRESS(192,168,1,33)

	// read/write int variable
	unsigned int xdriver_getValue(const char *tagName);
	bool xdriver_setValue(const char *tagName, const unsigned int value);

	// read/write float variable
	float xdriver_getValue_float(const char *tagName);
	bool xdriver_setValue_float(const char *tagName, const float value);

	// read/write in memory format
	bool xdriver_get(const char *tagName,void *pData,unsigned int len);
	bool xdriver_set(const char *tagName,const void *pData,unsigned int len);

	// donot use, will delete next version
	unsigned int xdriver_struct_get(const char *tagName,void *structData);
	bool xdriver_struct_set(const char *tagName, void *structData);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
