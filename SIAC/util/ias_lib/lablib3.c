/*========================================================================*/
/*                                                                        */
/*                         PDS Label Library Lite                         */
/*                               (lablib3)                                */
/*                                                                        */
/*  Version:                                                              */
/*                                                                        */
/*      1.0Beta    Mar 31, 1994                                           */
/*      1.0        Jan 23, 1995                                           */
/*      1.1        Feb 23, 1995                                           */
/*      1.2        Jun 06, 1995  Preliminary                              */
/*                                                                        */
/*  Change History:                                                       */
/*                                                                        */
/*      03-31-94    Original code                                         */
/*	01-09-95    jsh - Changed OBJECT to OBJDESC                       */
/*	01-09-95    jsh - Corrected strcmp and strncmp == NULL            */
/*	02-16-95    jsh - Applied LASP changes (from OA dvlp - s. monk)   */
/*			- Function Prototypes                             */
/*			- Filename Units                                  */
/*			- Filename for SUN/Unix                           */
/*			- TB_MAX_BUFFER                                   */
/*			- Several 0 -> NULL in function calls             */
/*	02-20-95    jsh - Added OdlPrintLine (reduced # fprintf)          */
/*	06-06-95    jsh/gmw - Allow SFDU without "= SFDU"                 */
/*	06-06-95    jsh/gmw - Stop gap for '/' followed by '*' in text strings */
/*                                                                        */
/*                                                                        */
/*========================================================================*/

#include "lablib3.h"

long odl_message_count = {0};





/*========================================================================*/
/*                                                                        */
/*                          Label Parse routines                          */
/*                                                                        */
/*========================================================================*/



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlParseLabelFile                                               */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine causes a file containing a PDS label to be parsed  */
/*      and its "^" keywords expanded.  It returns a pointer to the     */
/*      root of the OBJECT tree data structure.                         */
/*                                                                      */
/*      The "filespec" parameter is the full path and name of the       */
/*      label file to be parsed.                                        */
/*                                                                      */
/*      The "message_fname" parameter is the name of the file where     */
/*      parser error messages will be written.  If the file exists,     */
/*      the messages are appended, otherwise a new file is created.     */
/*      A NULL value passed in causes messages to be sent to stdout.    */
/*                                                                      */
/*      The "suppress_messages" parameter is a flag that tells the code */
/*      whether or not to print parser error messages.  A value of TRUE */
/*      (1) tells the code to supress all messages.  If this parameter  */
/*      is 1, it doesn't matter what you specified with the             */
/*      "message_fname".  Nothing will be written to that file.  If a   */
/*      zero (0) is passed in, then messages will be written to the     */
/*      file you specified with the "message_fname" parameter.          */
/*                                                                      */
/*      The expand parameter is a flag that controls whether or not     */
/*      ^STRUCTURE and ^CATALOG keywords are expanded.  To expand one   */
/*      of these keywords means to take the file name it points to,     */
/*      parse its contents, and insert the results into the tree right  */
/*      where the keyword is sitting.  In other words, these keywords   */
/*      function just like include files in "C".  These are the values  */
/*      that can be passed in:                                          */
/*                                                                      */
/*             ODL_EXPAND_STRUCTURE - expand ^STRUCTURE keywords only   */
/*             ODL_EXPAND_CATALOG   - expand ^CATALOG keywords only     */
/*             ODL_EXPAND_STRUCTURE | ODL_EXPAND_CATALOG - expand       */
/*                   both keywords (the "|" character is the logical    */
/*                   "or" of both values).                              */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*      WARNING:  The value returned by this routine points to memory   */
/*                allocated by this routine (sometimes quite a bit of   */
/*                memory!).  Be sure to deallocate it using the         */
/*                OdlFreeTree routine.                                  */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlParseLabelFile (filespec, message_fname, expand, suppress_messages)

char *filespec;
char *message_fname;
MASK expand;
int suppress_messages;

#else

OBJDESC *OdlParseLabelFile (char *filespec, char *message_fname, MASK expand, 
                            int suppress_messages)

#endif

{
    OBJDESC *root = {NULL};
    
    root = (OBJDESC *) OdlParseFile(filespec,NULL,message_fname,NULL,suppress_messages,1,1,0);
    root = (OBJDESC *) OdlExpandLabelFile(root, message_fname, expand,
                                         suppress_messages);

    return(root);

}  /*  End:  "OdlParseLabelFile"  */




/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlParseLabelString                                             */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine causes a character string containing an ODL        */
/*      statement to be parsed and its "^" keywords expanded.  It       */
/*      returns a pointer to the root of the OBJECT tree data structure.*/
/*                                                                      */
/*      WARNING:  The value returned by this routine points to memory   */
/*                allocated by this routine.  Be sure to deallocate it  */
/*                using the OdlFreeTree routine.                        */
/*                                                                      */
/*      WARNING:  This routine will try to create a temporary file in   */
/*                the following locations, depending on the system:     */
/*                                                                      */
/*                UNIX:        ~/<tmp fname>.tmp                        */
/*                VMS:         sys$login:<tmp fname>.tmp                */
/*                MSDOS:       C:\<tmp fname>.tmp                       */
/*                All others:  <tmp fname>.tmp                          */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlParseLabelString (odl_string, message_fname, 
                             expand, suppress_messages)

char *odl_string;
char *message_fname;
MASK expand;
int suppress_messages;

#else

OBJDESC *OdlParseLabelString (char *odl_string, char *message_fname, 
                             MASK expand, int suppress_messages)

#endif
{
    OBJDESC *root = {NULL};
    FILE *tmp_fptr = {NULL};
    char *tmp_fname = {NULL};

    tmp_fname = (char *) OdlTempFname();

    if (tmp_fname == NULL)
    {
        OdlPrintMessage(message_fname,NULL,0,
                        "Unable to create a temporary file", suppress_messages);
    }

    if ((tmp_fptr = (FILE *) fopen(tmp_fname, "w")) != NULL)
    {
        (void)fprintf(tmp_fptr, "%s", odl_string);
        (void)fclose(tmp_fptr);
        root = (OBJDESC *) OdlParseLabelFile(tmp_fname, message_fname, 
                                            expand, suppress_messages);
#ifdef VMS
        AppendString(tmp_fname, ";*")
#endif
        (void)remove(tmp_fname);
    }

    LemmeGo(tmp_fname)
    
    return(root);

}  /*  End:  "OdlParseLabelString"  */




/*========================================================================*/
/*                                                                        */
/*                        Label Expand routines                           */
/*                                                                        */
/*========================================================================*/

/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlExpandLabelFile                                              */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates "^STRUCTURE" and "^CATALOG" keywords and   */
/*      expands them.  To "expand" means to extract the file name       */
/*      pointed to by the keyword, parse its contents, and insert the   */
/*      resulting tree into the label right where the keyword is        */
/*      sitting.                                                        */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlExpandLabelFile (object, message_fname, expand, suppress_messages)

OBJDESC *object;
char *message_fname;
MASK expand;
int suppress_messages;

#else

OBJDESC *OdlExpandLabelFile (OBJDESC *object, char *message_fname, MASK expand, 
                            int suppress_messages)

#endif
{
    KEYWORD *kwd = {NULL};
    KEYWORD *new_kwd = {NULL};
    KEYWORD *save_kwd = {NULL};
    OBJDESC *temp_root = {NULL};
    OBJDESC *new_obj = {NULL};
    OBJDESC *save_obj = {NULL};
    FILE *l_ptr = {NULL};
    unsigned long start_loc = {1};
    unsigned short loc_type = {ODL_RECORD_LOCATION};
    unsigned short done = {FALSE};
    char *fspec = {NULL};
    char *fname = {NULL};
    char *keyword_name = {NULL};
    char error_message[5*(TB_MAXLINE + TB_MAXPATH + TB_MAXFNAME)];

    /*  Let's expand all ^STRUCTURE keywords, shall we?  */
    if ((expand&ODL_EXPAND_STRUCTURE) == ODL_EXPAND_STRUCTURE)
    {
        expand -= ODL_EXPAND_STRUCTURE;
        CopyString(keyword_name, "^*STRUCTURE")
    }
    else
    /*  On second thought, let's expand all ^CATALOG keywords  */
    if ((expand&ODL_EXPAND_CATALOG) == ODL_EXPAND_CATALOG)
    {
        expand -= ODL_EXPAND_CATALOG;
        CopyString(keyword_name, "^*CATALOG")
    }
    else
    /*  Hmmm. I guess we have nothing left to expand.  */
    {
        expand = ODL_NOEXPAND;
        done = TRUE;
    }

    /*  Keep expanding until we can expand no more forever  */
    while (! done)
    {
        /*  Find the expand keyword wherever it my be hiding  */
        kwd = (KEYWORD *) OdlFindKwd(object, keyword_name,
                                     NULL, 1, ODL_RECURSIVE_DOWN);

        /*  We're done if there aren't any more keywords to expand  */
        if (kwd == NULL)
            done = TRUE;
        else
        {
            /*  Get the file name, minus quotes and blanks, sans path  */
            fname = (char *) OdlGetFileName(kwd, &start_loc, &loc_type);

            /*  We're in trouble if we've encountered this file before  */
            if (ExpandIsRecursive(kwd, fname))
            {
                (void)sprintf(error_message, 
                        "Recursive %s statement found in file:  %s", 
                        keyword_name, fname);
                OdlPrintMessage(message_fname,NULL,(long)kwd->line_number,
                                error_message, suppress_messages);
            }
            else
            {
                /*  Figure out exactly where the file is located  */
                fspec = (char *) OdlGetFileSpec(fname);

                /*  We're in trouble if we can't find the file  */
                if (fspec == NULL)
                {
                    (void)sprintf(error_message, 
                            "Unable to locate %s file:  %s",
                            keyword_name, fname);
                    OdlPrintMessage(message_fname,NULL,(long)kwd->line_number,
                                    error_message, suppress_messages);
                }
                else
                {
                    l_ptr = (FILE *) OdlLocateStart(fspec, start_loc, loc_type);

   
                    /*  Parse the file  */
                    temp_root = (OBJDESC *) OdlParseFile(fspec,l_ptr,
                                                        message_fname,NULL,0,1,1,1);
                    
                    /*  Was there anything in the file to parse?  */
                    if (temp_root != NULL)
                    {
                        /*  Append any keywords  */
                        for (new_kwd=temp_root->first_keyword;
                                 new_kwd != NULL; new_kwd = save_kwd)
                        {
                            save_kwd = new_kwd->right_sibling;
                            (void)OdlPasteKwd((KEYWORD *) OdlCutKwd(new_kwd),
                                        kwd->parent);
                        }

                        /*  Append any sub-objects  */
                        for (new_obj=temp_root->first_child;
                                 new_obj != NULL; new_obj = save_obj)
                        {
                            save_obj = new_obj->right_sibling;
                            (void)OdlPasteObjDesc((OBJDESC *) OdlCutObjDesc(new_obj),
                                            kwd->parent);
                        }

                        /*  Deallocate the temporary root  */
                        temp_root->first_keyword = NULL;
                        temp_root->first_child = NULL;
                        temp_root = (OBJDESC *) OdlFreeTree(temp_root);

                    }  /*  End:  "if (temp_root != NULL) ..."  */
    
                    /*  Free the file spec storage  */
                    LemmeGo(fspec)

                }  /*  End:  "if (fspec == NULL) ... else ..."  */

            }  /*  End:  "if (ExpandIsRecursive( ... else ..."  */

            (void)OdlFreeKwd((KEYWORD *)OdlCutKwd(kwd));

            /*  Free the file name storage  */
            LemmeGo(fname)
            
        }  /*  End:  "if (kwd == NULL) ... else ..."  */

    }  /*  End:  "while (! done) ..."  */

    /*  Free the keyword name storage  */
    LemmeGo(keyword_name)

    /*  Check and see if there are any other keywords to expand  */
    if (expand != ODL_NOEXPAND)
    {
        object = (OBJDESC *) OdlExpandLabelFile(object, message_fname, 
                                               expand, suppress_messages);
    }

    /*  Return the root of the expanded tree  */
    return(object);

}  /*  End:  "OdlExpandLabelFile"  */



/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

unsigned short ExpandIsRecursive (keyword, exp_fname)

KEYWORD *keyword;
char *exp_fname;

#else

unsigned short ExpandIsRecursive (KEYWORD *keyword, char *exp_fname)

#endif
{
    OBJDESC *obj = {NULL};
    char *temp_fname = {NULL};
    unsigned short found = {FALSE};

    if ((keyword != NULL) && (exp_fname != NULL))
    {
#if (defined( VAX) || defined( ALPHA_VMS))
        UpperCase(exp_fname)
#endif

        CopyString(temp_fname, keyword->file_name)

#if (defined( VAX) || defined( ALPHA_VMS))
        UpperCase(temp_fname)
#endif

        found = (strcmp(temp_fname, exp_fname) == 0);
        LemmeGo(temp_fname)

        for (obj=keyword->parent; 
                  ((! found) && (obj != NULL)); obj=obj->parent)
        {
            CopyString(temp_fname, obj->file_name)

#if (defined( VAX) || defined( ALPHA_VMS))
            UpperCase(temp_fname)
#endif

            found = (strcmp(temp_fname, exp_fname) == 0);
            LemmeGo(temp_fname)
        }
    }

    return(found);

}  /*  End:  "ExpandIsRecursive"  */





/*========================================================================*/
/*                                                                        */
/*                     Object description routines                        */
/*                                                                        */
/*========================================================================*/

/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlFindObjDesc                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates an object within a parsed label by its     */
/*      class name (like TABLE), by its position (look for the seventh  */
/*      table object in the label), by a particular keyword present     */
/*      in the object (like NAME), or by a particular value that a      */
/*      particular keyword has (like START_BYTE = 76).                  */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlFindObjDesc(start_object, object_class, keyword_name, 
                       keyword_value, object_position, search_scope)

OBJDESC *start_object;
const char *object_class;
const char *keyword_name;
char *keyword_value;
unsigned long object_position;
unsigned short search_scope;

#else

OBJDESC *OdlFindObjDesc( OBJDESC *start_object, const char *object_class, 
                        const char *keyword_name, char *keyword_value, 
                        unsigned long object_position, 
                        unsigned short search_scope)

#endif
{
    OBJDESC *found_object = {NULL};
    OBJDESC *obj = {NULL};
    KEYWORD *kwd = {NULL};
    unsigned short found = {FALSE};
    unsigned short scope = {search_scope};
    unsigned long current_position = {0};

    for (obj=start_object;
          ((obj != NULL) && (! found));
            obj = (OBJDESC *) OdlNextObjDesc(obj, start_object->level, &scope))
    {
        if (object_class == NULL)
             found = TRUE;
        else
             found = OdlWildCardCompare(object_class, obj->class);

        if ((found) && (keyword_name != NULL))
        {
            kwd = (KEYWORD *) OdlFindKwd(obj, keyword_name, 
                                         NULL, 1, ODL_THIS_OBJECT);
            found = (kwd != NULL);
        }

        if ((found) && (keyword_value != NULL))
            found = OdlWildCardCompare(keyword_value, (char *) OdlGetKwdValue(kwd));

        if ((found) && (object_position != 0))
            found = ((++current_position) == object_position);

        if (found) found_object = obj;    

    }  /*  End:  "for (obj=start_object; ..."  */

    return(found_object);

}  /*  End:  "OdlFindObjDesc"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlNextObjDesc                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the next object in the tree based on the   */
/*      search_scope passed in.                                         */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlNextObjDesc (object, root_level, search_scope)

OBJDESC *object;
unsigned long root_level;
unsigned short *search_scope;

#else

OBJDESC *OdlNextObjDesc (OBJDESC *object, unsigned long root_level, 
                        unsigned short *search_scope)

#endif
{
    OBJDESC *next_object = {NULL};

    if (object != NULL)
    {
        switch (*search_scope)
        {
            /*  look only in the current object  */
            case ODL_THIS_OBJECT    :  next_object = NULL;
                                       break;

            /*  look at the current object's first child now, and its   */
            /*  child's right siblings in subsequent searches           */
            case ODL_CHILDREN_ONLY  :  next_object = object->first_child;
                                       *search_scope = ODL_SIBLINGS_ONLY;
                                       break;

            /*  look at the current object's right sibling  */
            case ODL_SIBLINGS_ONLY  :  next_object = object->right_sibling;
                                       break;

            /*  treat the current object as the root of a sub-tree  */
            case ODL_RECURSIVE_DOWN :  next_object = (OBJDESC *) OdlTraverseTree(object, root_level);
                                       break;

            /*  search children, then siblings, then move up to parent  */
            /*  keep going until the end of the label is reached        */
            default                 :  next_object = (OBJDESC *) 
                                                      OdlTraverseTree(object, 
                                                                     (unsigned long) 0);
                                       break;

        }  /*  End:  "switch (*search_scope) ..."  */

    }  /*  End:  "if (object != NULL) ..."  */

    return(next_object);

}  /*  End:  "OdlNextObjDesc"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlCutObjDesc                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine cuts an ODL object structure out of an ODL tree    */
/*      and returns a pointer to it.  All references to it in the tree  */
/*      are removed, and all references to the original tree within the */
/*      object are removed.                                             */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlCutObjDesc (object)

OBJDESC *object;

#else

OBJDESC *OdlCutObjDesc (OBJDESC *object)

#endif
{
    if (object != NULL)
    {
        if (object->right_sibling == NULL)
            object->parent->last_child = object->left_sibling;
        else
            object->right_sibling->left_sibling = object->left_sibling;

        if (object->left_sibling == NULL)
            object->parent->first_child = object->right_sibling;
        else
            object->left_sibling->right_sibling = object->right_sibling;

        object->parent = NULL;
        object->left_sibling = NULL;
        object->right_sibling = NULL;

    }  /*  End:  "if (object != NULL) ..."  */
    
    return(object);

}  /*  End routine:  "OdlCutObjDesc"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPasteObjDesc                                                 */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine adds an object to a tree as the last child of the  */
/*      parent_object.                                                  */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlPasteObjDesc (new_object, parent_object)

OBJDESC *new_object;
OBJDESC *parent_object;

#else

OBJDESC *OdlPasteObjDesc (OBJDESC *new_object, OBJDESC *parent_object)

#endif
{
    if ((new_object != NULL) && (parent_object != NULL))
    {
        new_object->left_sibling = parent_object->last_child;
        new_object->right_sibling = NULL;
        new_object->parent = parent_object;
    
        if (parent_object->first_child == NULL)
            parent_object->first_child = new_object;

        if (parent_object->last_child != NULL)
            parent_object->last_child->right_sibling = new_object;

        parent_object->last_child = new_object;

        OdlAdjustObjDescLevel(new_object);

    }  /*  End:  "if ((new_object != NULL) && ..."  */
    
    return(new_object);

}  /*  End routine:  "OdlPasteObjDesc"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPasteObjDescBefore                                           */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine adds an object to a tree as the left sibling of    */
/*      the old_object.                                                 */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlPasteObjDescBefore (new_object, old_object)

OBJDESC *new_object;
OBJDESC *old_object;

#else

OBJDESC *OdlPasteObjDescBefore (OBJDESC *new_object, OBJDESC *old_object)

#endif
{
    if ((new_object != NULL) && (old_object != NULL))
    {
        new_object->left_sibling = old_object->left_sibling;
        new_object->right_sibling = old_object;
        new_object->parent = old_object->parent;
    
        if (old_object->left_sibling == NULL)
            old_object->parent->first_child = new_object;
        else
            old_object->left_sibling->right_sibling = new_object;

        old_object->left_sibling = new_object;

        OdlAdjustObjDescLevel(new_object);

    }  /*  End:  "if ((new_object != NULL) && ..."  */
    
    return(new_object);

}  /*  End routine:  "OdlPasteObjDescBefore"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPasteObjDescAfter                                            */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine adds an object to a tree as the right sibling of   */
/*      the old_object.                                                 */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlPasteObjDescAfter (new_object, old_object)

OBJDESC *new_object;
OBJDESC *old_object;

#else

OBJDESC *OdlPasteObjDescAfter (OBJDESC *new_object, OBJDESC *old_object)

#endif
{
    if ((new_object != NULL) && (old_object != NULL))
    {
        new_object->right_sibling = old_object->right_sibling;
        new_object->left_sibling = old_object;
        new_object->parent = old_object->parent;
    
        if (old_object->right_sibling == NULL)
            old_object->parent->last_child = new_object;
        else
            old_object->right_sibling->left_sibling = new_object;

        old_object->right_sibling = new_object;

        OdlAdjustObjDescLevel(new_object);

    }  /*  End:  "if ((new_object != NULL) && ..."  */
    
    return(new_object);

}  /*  End routine:  "OdlPasteObjDescAfter"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlCopyObjDesc                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine makes a copy of an object and returns a pointer    */
/*      to the copy.  All fields are duplicated except for references   */
/*      to the original tree, which are removed.                        */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlCopyObjDesc (object)

OBJDESC *object;

#else

OBJDESC *OdlCopyObjDesc (OBJDESC *object)

#endif
{
    OBJDESC *new_object = {NULL};

    if (object != NULL)
    {
        new_object = OdlNewObjDesc(object->class, 
                               object->pre_comment, object->line_comment,
                               object->post_comment, object->end_comment,
                               object->file_name, object->is_a_group, 
                               (long)object->line_number);
    }
    
    return(new_object);

}  /*  End routine:  "OdlCopyObjDesc"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlNewObjDesc                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine creates a new object structure and initializes     */
/*      its fields with the values passed in.                           */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlNewObjDesc (object_class, pre_comment, line_comment, post_comment, 
                       end_comment, file_name, is_a_group, line_number)
const char *object_class;
const char *pre_comment;
const char *line_comment;
const char *post_comment;
const char *end_comment;
const char *file_name;
short is_a_group;
long line_number;

#else

OBJDESC *OdlNewObjDesc (const char *object_class, const char *pre_comment, 
                       const char *line_comment, const char *post_comment, 
                       const char *end_comment, const char *file_name, 
                       short is_a_group, long line_number)

#endif
{
    OBJDESC *new_object = {NULL};

    if ((new_object = (OBJDESC *)malloc(sizeof(OBJDESC))) == NULL)
        SayGoodbye()
    else
    {
        CopyString(new_object->class, object_class)
        CopyString(new_object->pre_comment, pre_comment)
        CopyString(new_object->line_comment, line_comment)
        CopyString(new_object->post_comment, post_comment)
        CopyString(new_object->end_comment, end_comment)
        CopyString(new_object->file_name, file_name)

        new_object->is_a_group = is_a_group;
        new_object->child_count = 0;
        new_object->line_number = line_number;
        new_object->level = 0;
        new_object->parent = NULL;
        new_object->left_sibling = NULL;
        new_object->right_sibling = NULL;
        new_object->first_child = NULL;
        new_object->last_child = NULL;
        new_object->first_keyword = NULL;
        new_object->last_keyword = NULL;
        new_object->appl1 = NULL;
        new_object->appl2 = NULL;

    }  /*  End:  "if ((new_object = ... else ..."  */

    return(new_object);

}  /*  End routine:  "OdlNewObjDesc"  */




/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetLabelVersion                                              */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a pointer to a character string containing */
/*      the ODL version of the label.  It looks for this information in */
/*      the ODL_VERSION_NUMBER keyword.                                 */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO THE ACTUAL VALUE OF THE         */
/*                ODL_VERSION_NUMBER KEYWORD AND MUST NOT BE FREED.     */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlGetLabelVersion (object)

OBJDESC *object;

#else

char *OdlGetLabelVersion (OBJDESC *object)

#endif
{
    KEYWORD *kwd = {NULL};
    char *version = {NULL};

    if (object != NULL)
    {
        kwd = (KEYWORD *) OdlFindKwd(object, "PDS_VERSION_ID", 
                                     NULL, 1, ODL_THIS_OBJECT);
        if (kwd == NULL) 
        {
            kwd = (KEYWORD *) OdlFindKwd(object, "ODL_VERSION_NUMBER", 
                                         NULL, 1, ODL_THIS_OBJECT);
	}

        if (kwd != NULL) 
            version = (char *) OdlGetKwdValue(kwd);
    }

    return(version);

}  /*  End:  "OdlGetLabelVersion"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetObjDescClassName                                          */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns the class name of an object.               */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO THE ACTUAL INFORMATION STORED   */
/*                IN THE ODL OBJECT STRUCTURE AND MUST NOT BE FREED.    */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlGetObjDescClassName (object)

OBJDESC *object;

#else

char *OdlGetObjDescClassName (OBJDESC *object)

#endif
{
    char *class_name = {NULL};

    if (object != NULL)
        class_name = object->class;

    return(class_name);

}  /*  End:  "OdlGetObjDescClassName"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetObjDescChildCount                                         */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a count of the immediate children of an    */
/*      object.  It does not count children of children, etc.           */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

int OdlGetObjDescChildCount (object)

OBJDESC *object;

#else

int OdlGetObjDescChildCount (OBJDESC *object)

#endif
{
    OBJDESC *obj = {NULL};
    int child_count = {0};

    if (object != NULL)
    {
        for (obj=object->first_child; obj != NULL; obj=obj->right_sibling)
            ++child_count;
    }

    return(child_count);

}  /*  End:  "OdlGetObjDescChildCount"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetObjDescLevel                                              */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns the nesting level of an object.  The ROOT  */
/*      object in a tree is always defined to be level 0.               */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

int OdlGetObjDescLevel (object)

OBJDESC *object;

#else

int OdlGetObjDescLevel (OBJDESC *object)

#endif
{
    OBJDESC *obj = {NULL};
    int level = {0};

    if (object != NULL)
    {
        for (obj=object->parent; obj != NULL; obj=obj->parent)
            ++level;
    }

    return(level);

}  /*  End:  "OdlGetObjDescLevel"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlAdjustObjDescLevel                                           */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine changes the nesting level of an object and all of  */
/*      its subobjects so they fit in with their place in the overall   */
/*      ODL tree.  This is particularly useful when objects are cut     */
/*      from one tree and pasted into another tree, perhaps higher or   */
/*      lower in the nesting hierarchy then they were in the original   */
/*      tree.                                                           */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

void OdlAdjustObjDescLevel (object)

OBJDESC *object;

#else

void OdlAdjustObjDescLevel (OBJDESC *object)

#endif
{
    OBJDESC *obj = {NULL};
    unsigned short scope = {ODL_RECURSIVE_DOWN};

    for (obj=object; obj != NULL; 
             obj = (OBJDESC *) OdlNextObjDesc(obj, object->level, &scope))
    {
        obj->level = (obj->parent == NULL) ? 0 : (1 + obj->parent->level);
    }

    return;

}  /*  End routine:  "OdlAdjustObjDescLevel"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetObjDescParent                                             */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a pointer to an object's parent.           */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO AN EXISTING ODL OBJECT AND      */
/*                AND MUST NOT BE FREED.                                */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlGetObjDescParent (object)

OBJDESC *object;

#else

OBJDESC *OdlGetObjDescParent (OBJDESC *object)

#endif
{
    OBJDESC *parent = {NULL};

    if (object != NULL)
        parent = object->parent;

    return(parent);

}  /*  End:  "OdlGetObjDescParent"  */




/*========================================================================*/
/*                                                                        */
/*                           Keyword routines                             */
/*                                                                        */
/*========================================================================*/


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlFindKwd                                                      */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the keyword in a label that satisfies the  */
/*      requirements passed in:  The object where the search is to      */
/*      begin, the name of the keyword, a particular value that the     */
/*      keyword must have, which version of the keyword we want (if     */
/*      there are duplicates), and the search scope we want to use to   */
/*      limit the objects searched.                                     */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlFindKwd (start_object, keyword_name, keyword_value, 
                         keyword_position, search_scope)
OBJDESC *start_object;
const char *keyword_name;
char *keyword_value;
unsigned long keyword_position;
unsigned short search_scope;

#else

KEYWORD *OdlFindKwd (OBJDESC *start_object, const char *keyword_name, 
                     char *keyword_value, unsigned long keyword_position, 
                     unsigned short search_scope)

#endif
{
    OBJDESC *obj = {NULL};
    KEYWORD *kwd = {NULL};
    KEYWORD *found_kwd = {NULL};
    unsigned short found = {FALSE};
    unsigned short scope = {search_scope};
    unsigned long current_position = {0};

    for (obj=start_object;
            ((obj != NULL) && (! found));
                obj = (OBJDESC *) OdlNextObjDesc(obj, start_object->level, &scope))
    {
        for (kwd=obj->first_keyword; ((kwd != NULL) && (! found)); kwd=kwd->right_sibling)
        {
            if (keyword_name == NULL)
                found = TRUE;
            else
                found = OdlWildCardCompare(keyword_name, kwd->name);

            if ((found) && (keyword_value != NULL))
                found = OdlWildCardCompare(keyword_value, (char *) OdlGetKwdValue(kwd));
    
            if ((found) && (keyword_position != 0))
                found = ((++current_position) == keyword_position);
    
            if (found) found_kwd = kwd;    

        }  /*  End:  "for (kwd=obj-> ..."  */

    }  /*  End:  "for (obj=start_object; ..."  */

    return(found_kwd);

}  /*  End:  "OdlFindKwd"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlNextKwd                                                      */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the keyword in a label that satisfies the  */
/*      requirements passed in:  The object where the search is to      */
/*      begin, the name of the keyword, a particular value that the     */
/*      keyword must have, which version of the keyword we want (if     */
/*      there are duplicates), and the search scope we want to use to   */
/*      limit the objects searched.                                     */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlNextKwd (start_keyword, keyword_name, keyword_value, 
                     keyword_position, search_scope)

KEYWORD *start_keyword;
char *keyword_name;
char *keyword_value;
unsigned long keyword_position;
unsigned short search_scope;

#else

KEYWORD *OdlNextKwd (KEYWORD *start_keyword, char *keyword_name, 
                     char *keyword_value, unsigned long keyword_position, 
                     unsigned short search_scope)

#endif
{
    OBJDESC *start_object = {NULL};
    OBJDESC *obj = {NULL};
    KEYWORD *kwd = {NULL};
    KEYWORD *found_kwd = {NULL};
    unsigned short found = {FALSE};
    unsigned short scope = {search_scope};
    unsigned long current_position = {0};

    if (start_keyword != NULL)
    {
        start_object = start_keyword->parent;
        obj = start_object;
        kwd = start_keyword; 
    
        do
        {
            for ( ; ((kwd != NULL) && (! found)); kwd=kwd->right_sibling)
            {
                if (keyword_name == NULL)
                    found = TRUE;
                else
                    found = OdlWildCardCompare(keyword_name, kwd->name);
        
                if ((found) && (keyword_value != NULL))
                    found = OdlWildCardCompare(keyword_value, (char *) OdlGetKwdValue(kwd));
        
                if ((found) && (keyword_position != 0))
                    found = ((++current_position) == keyword_position);
        
                if (found) found_kwd = kwd;    
        
            }  /*  End:  "for (kwd=start_keyword; ..."  */
    
            if (! found)
            {
                obj = (OBJDESC *) OdlNextObjDesc(obj, start_object->level, &scope);
                kwd = (KEYWORD *) OdlGetFirstKwd(obj);
            }
    
        }  while ((obj != NULL) && (! found));

    }  /*  End:  "if (start_keyword != NULL) ..."  */

    return(found_kwd);

}  /*  End:  "OdlNextKwd"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlCutKwd                                                       */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine removes a keyword from an object and returns a     */
/*      pointer to it.  All references to the object within the keyword */
/*      are removed, and all references to the keyword within the       */
/*      object are removed.                                             */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO AN EXISTING KEYWORD STRUCTURE   */
/*                AND MUST NOT BE FREED.                                */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlCutKwd (keyword)

KEYWORD *keyword;

#else

KEYWORD *OdlCutKwd (KEYWORD *keyword)

#endif
{
    if (keyword != NULL)
    {
        if (keyword->left_sibling != NULL)
            keyword->left_sibling->right_sibling = keyword->right_sibling;

        if (keyword->right_sibling != NULL)
            keyword->right_sibling->left_sibling = keyword->left_sibling;

        if (keyword->parent->first_keyword == keyword)
            keyword->parent->first_keyword = keyword->right_sibling;

        if (keyword->parent->last_keyword == keyword)
            keyword->parent->last_keyword = keyword->left_sibling;

        keyword->parent = NULL;
        keyword->left_sibling = NULL;
        keyword->right_sibling = NULL;

    }  /*  End:  "if ((keyword != NULL) && ..."  */

    return(keyword);

}  /*  End routine:  "OdlCutKwd"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPasteKwd                                                     */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine adds a keyword to the end of an object's keyword   */
/*      list.                                                           */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlPasteKwd (keyword, object)

KEYWORD *keyword;
OBJDESC *object;

#else

KEYWORD *OdlPasteKwd (KEYWORD *keyword, OBJDESC *object)

#endif
{
    if ((keyword != NULL) && (object != NULL))
    {
        keyword->parent = object;
        keyword->left_sibling = object->last_keyword;
        keyword->right_sibling = NULL;

        if (object->first_keyword == NULL)
            object->first_keyword = keyword;

        if (object->last_keyword != NULL)
            object->last_keyword->right_sibling = keyword;

        object->last_keyword = keyword;

    }  /*  End:  "if ((keyword != NULL) && ..."  */

    return(keyword);

}  /*  End routine:  "OdlPasteKwd"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPasteKwdBefore                                               */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine adds a keyword to an object as the left sibling of */
/*      the old_keyword.                                                */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlPasteKwdBefore (new_keyword, old_keyword)

KEYWORD *new_keyword;
KEYWORD *old_keyword;

#else

KEYWORD *OdlPasteKwdBefore (KEYWORD *new_keyword, KEYWORD *old_keyword)

#endif
{
    if ((new_keyword != NULL) && (old_keyword != NULL))
    {
        new_keyword->parent = old_keyword->parent;
        new_keyword->left_sibling = old_keyword->left_sibling;
        new_keyword->right_sibling = old_keyword;

        if (old_keyword->left_sibling == NULL)
            old_keyword->parent->first_keyword = new_keyword;
        else
            old_keyword->left_sibling->right_sibling = new_keyword;

        old_keyword->left_sibling = new_keyword;

    }  /*  End:  "if ((new_keyword != NULL) && ..."  */

    return(new_keyword);

}  /*  End routine:  "OdlPasteKwdBefore"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPasteKwdAfter                                                */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine adds a keyword to an object as the right sibling   */
/*      of the old_keyword.                                             */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlPasteKwdAfter (new_keyword, old_keyword)

KEYWORD *new_keyword;
KEYWORD *old_keyword;

#else

KEYWORD *OdlPasteKwdAfter (KEYWORD *new_keyword, KEYWORD *old_keyword)

#endif
{
    if ((new_keyword != NULL) && (old_keyword != NULL))
    {
        new_keyword->parent = old_keyword->parent;
        new_keyword->right_sibling = old_keyword->right_sibling;
        new_keyword->left_sibling = old_keyword;

        if (old_keyword->right_sibling == NULL)
            old_keyword->parent->last_keyword = new_keyword;
        else
            old_keyword->right_sibling->left_sibling = new_keyword;

        old_keyword->right_sibling = new_keyword;

    }  /*  End:  "if ((new_keyword != NULL) && ..."  */

    return(new_keyword);

}  /*  End routine:  "OdlPasteKwdAfter"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlCopyKwd                                                      */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine makes a copy of a keyword and returns a pointer to */
/*      it.  All of the keyword's fields are duplicated except for      */
/*      references to the parent object, which are removed.             */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlCopyKwd (keyword)

KEYWORD *keyword;

#else

KEYWORD *OdlCopyKwd (KEYWORD *keyword)

#endif
{
    KEYWORD *new_keyword = {NULL};

    if (keyword != NULL)
    {
        new_keyword = OdlNewKwd(keyword->name, keyword->value, 
                                 keyword->pre_comment, keyword->line_comment,
                                 keyword->file_name, (long)keyword->line_number);
    }

    return(new_keyword);

}  /*  End routine:  "OdlCopyKwd"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlNewKwd                                                       */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine creates a new keyword structure, initializes       */
/*      its fields with the values passed in, and returns a pointer     */
/*      to it.                                                          */
/*                                                                      */	
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlNewKwd (keyword_name, value_text, pre_comment, 
                    line_comment, file_name, line_number)

char *keyword_name;
char *value_text;
char *pre_comment;
char *line_comment;
char *file_name;
long line_number;

#else

KEYWORD *OdlNewKwd (char *keyword_name, char *value_text, char *pre_comment, 
                    char *line_comment, char *file_name, long line_number)

#endif
{
    KEYWORD *new_keyword = {NULL};

    if ((new_keyword = (KEYWORD *)malloc(sizeof(KEYWORD))) == NULL)
        SayGoodbye()
    else
    {
        CopyString(new_keyword->name, keyword_name)
        CopyString(new_keyword->pre_comment, pre_comment)
        CopyString(new_keyword->line_comment, line_comment)
        CopyString(new_keyword->file_name, file_name)
        CopyString(new_keyword->value, value_text)

        new_keyword->is_a_pointer = (keyword_name == NULL) ? FALSE : (*keyword_name == '^');

        if (value_text == NULL) 
        {
            new_keyword->size = 0;
            new_keyword->is_a_list = FALSE;
        }
        else
        {
            new_keyword->size = strlen(new_keyword->value);
            new_keyword->is_a_list = ((*value_text == '{') || (*value_text == '('));
        }

        new_keyword->line_number = line_number;
        new_keyword->parent = NULL;
        new_keyword->left_sibling = NULL;
        new_keyword->right_sibling = NULL;
        new_keyword->appl1 = NULL;
        new_keyword->appl2 = NULL;

    }  /*  End:  "if ((new_keyword = ... else ..."  */

    return(new_keyword);

}  /*  End routine:  "OdlNewKwd"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetFirstKwd                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a pointer to the first keyword data        */
/*      structure in an object definition structure.                    */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO THE ACTUAL INFORMATION STORED   */
/*                IN THE OBJECT DATA STRUCTURE AND MUST NOT BE FREED.   */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlGetFirstKwd (object)

OBJDESC *object;

#else

KEYWORD *OdlGetFirstKwd (OBJDESC *object)

#endif
{
    KEYWORD *kwd = {NULL};

    if (object != NULL)
        kwd = object->first_keyword;

    return(kwd);

}  /*  End:  "OdlGetFirstKwd"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetNextKwd                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a pointer to the next keyword data         */
/*      structure in an object definition's list of keyword structures. */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO THE ACTUAL INFORMATION STORED   */
/*                IN THE OBJECT DATA STRUCTURE AND MUST NOT BE FREED.   */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlGetNextKwd (keyword)

KEYWORD *keyword;

#else

KEYWORD *OdlGetNextKwd (KEYWORD *keyword)

#endif
{
    KEYWORD *kwd = {NULL};

    if (keyword != NULL)
        kwd = keyword->right_sibling;

    return(kwd);

}  /*  End:  "OdlGetNextKwd"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetKwdValue                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a pointer to a keyword's value.            */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO THE ACTUAL INFORMATION STORED   */
/*                IN THE KEYWORD DATA STRUCTURE AND MUST NOT BE FREED.  */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlGetKwdValue (keyword)

KEYWORD *keyword;

#else

char *OdlGetKwdValue (KEYWORD *keyword)

#endif
{
    char *value = {NULL};

    if (keyword != NULL)
        value = keyword->value;

    return(value);

}  /*  End:  "OdlGetKwdValue"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetAllKwdValues                                              */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine extracts each individual value from a set or       */
/*      sequence of values of a keyword, stores these values in a       */
/*      linked list, and returns a pointer to this list.                */
/*                                                                      */
/*      For example, if a keyword has this combination of sets and      */
/*      sequences as its value:                                         */
/*                                                                      */
/*         {red, (green, blue), {17, (("book.lbl", 345), orange)}}      */
/*                                                                      */
/*      Then the TB_STRING_LIST returned would contain:                 */
/*                                                                      */
/*         red                                                          */
/*         green                                                        */
/*         blue                                                         */
/*         17                                                           */
/*         "book.lbl"                                                   */
/*         345                                                          */
/*         orange                                                       */
/*                                                                      */
/*      WARNING:  The string list must be freed using the               */
/*                RemoveStringList macro (look in toolbox.h).           */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

TB_STRING_LIST *OdlGetAllKwdValues (keyword)

KEYWORD *keyword;

#else

TB_STRING_LIST *OdlGetAllKwdValues (KEYWORD *keyword)

#endif
{
    TB_STRING_LIST *value_list = {NULL};
    char *val_start = {NULL};
    char *val_stop = {NULL};
    char save_ch;

    if (keyword != NULL)
    {
        if (keyword->value != NULL)
        {
            for (val_start=(char *)OdlValueStart(keyword->value);
                     *val_start != '\0'; 
                         val_start=(char *)OdlValueStart(val_stop+1))
            {
                val_stop = (char *) OdlValueEnd(val_start);
                save_ch = *(val_stop + 1); *(val_stop + 1) = '\0';
                AddStringToList(val_start, value_list)
                *(val_stop + 1) = save_ch;
            }

        }  /*  End:  "if (keyword->value != NULL) ..."  */

    }  /*  End:  "if (keyword != NULL) ..."  */

    return(value_list);

}  /*  End:  "OdlGetAllKwdValues"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetKwdValueType                                              */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine determines the data type of a keyword's value and  */
/*      returns a data type symbolic id.  Possible symbolic ids are:    */
/*                                                                      */
/*         ODL_UNKNOWN   (can't tell what the heck it is)               */
/*         ODL_INTEGER   (handles optional leading plus or minus)       */
/*         ODL_REAL      (handles optional leading plus or minus,       */
/*                        scientific notation, and real exponents)      */
/*         ODL_SYMBOL    (unqouted or single quoted string of           */
/*                        characters)                                   */
/*         ODL_TEXT      (double quoted string of characters)           */
/*         ODL_DATE      (yyyy-mm-dd or yyyy-ddd)                       */
/*         ODL_DATE_TIME (yyyy-mm-ddThh:mm:ss.h)                        */
/*         ODL_SEQUENCE  (starts with a paren character "{")            */
/*         ODL_SET       (starts with a brace character "(")            */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

unsigned short OdlGetKwdValueType (keyword)

KEYWORD *keyword;

#else

unsigned short OdlGetKwdValueType (KEYWORD *keyword)

#endif
{
    unsigned short type = {ODL_UNKNOWN};
    if (keyword != NULL) type = (unsigned short) OdlDataType(keyword->value);
    return(type);

}  /*  End:  "OdlGetKwdValueType"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetKwdUnit                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the units part of a keyword's value,       */
/*      extracts it, stores it in a new character string, and returns   */
/*      a pointer to this new character string.                         */
/*                                                                      */
/*      WARNING:  This routine allocates memory for the return value    */
/*                that must be freed.                                   */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlGetKwdUnit (keyword)

KEYWORD *keyword;

#else

char *OdlGetKwdUnit (KEYWORD *keyword)

#endif
{
    char *c = {NULL};
    char *unit = {NULL};

    /*  If we were given a keyword to use  */
    if (keyword != NULL)
    {
        /*  Attempt to locate the units string  */
        c = (char *) strchr(keyword->value, '<');

        if (c != NULL)
        {
            /*  We found it!  Now copy it and make it upper case  */
            CopyString(unit, c)
            UpperCase(unit)
    
            /*  Close off the units string  */
            c = (char *) strchr(unit, '>');
            if (c != NULL) *(c + 1) = '\0';
        }

    }  /*  End:  "if (keyword != NULL) ..."  */

    return(unit);

}  /*  End:  "OdlGetKwdUnit"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetKwdName                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns the name of a keyword.                     */
/*                                                                      */
/*      WARNING:  NO MEMORY IS ALLOCATED BY THIS ROUTINE.  THE RETURN   */
/*                VALUE IS A POINTER TO THE ACTUAL INFORMATION STORED   */
/*                IN THE KEYWORD DATA STRUCTURE AND MUST NOT BE FREED.  */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlGetKwdName (keyword)

KEYWORD *keyword;

#else

char *OdlGetKwdName (KEYWORD *keyword)

#endif
{
    char *name = {NULL};

    if (keyword != NULL)
        name = keyword->name;

    return(name);

}  /*  End:  "OdlGetKwdName"  */





/*========================================================================*/
/*                                                                        */
/*                    Memory deallocation routines                        */
/*                                                                        */
/*========================================================================*/

/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlFreeTree                                                     */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine frees all memory used by an ODL tree. The return   */
/*      value is always NULL.                                           */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlFreeTree (object)

OBJDESC *object;

#else

OBJDESC *OdlFreeTree (OBJDESC *object)

#endif
{
    if (object != NULL)
    {
        (void)OdlFreeTree(object->first_child);
        (void)OdlFreeTree(object->right_sibling);
        (void)OdlFreeAllKwds(object);
        LemmeGo(object->class)      
        LemmeGo(object->pre_comment)
        LemmeGo(object->line_comment)
        LemmeGo(object->post_comment)
        LemmeGo(object->end_comment)
        LemmeGo(object->file_name)
        LemmeGo(object)
    }               
                    
    return(object);
                    
}  /*  End:  "OdlFreeTree"  */
                    
        


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlFreeAllKwds                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine frees all memory used by an object's keywords.     */
/*      When it's finished, all references to keywords are gone from    */
/*      the object.  The return value is always NULL.                   */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlFreeAllKwds (object)

OBJDESC *object;

#else

KEYWORD *OdlFreeAllKwds (OBJDESC *object)

#endif
{
    KEYWORD *kwd = {NULL};

    if (object != NULL)
    {
        for (kwd=object->first_keyword; kwd != NULL; 
                    kwd=(KEYWORD *) OdlFreeKwd(kwd)) ;

        object->first_keyword = NULL;
        object->last_keyword = NULL;
    }               
                    
    return(kwd);
                    
}  /*  End:  "OdlFreeAllKwds"  */
                    
        



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlFreeKwd                                                      */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine frees the memory used by a keyword.  The return    */
/*      value is always a pointer to the right sibling of the keyword.  */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

KEYWORD *OdlFreeKwd (keyword)

KEYWORD *keyword;

#else

KEYWORD *OdlFreeKwd (KEYWORD *keyword)

#endif
{
    KEYWORD *next_kwd = {NULL};

    if (keyword != NULL)
    {
        next_kwd = keyword->right_sibling;
        LemmeGo(keyword->name)        
        LemmeGo(keyword->file_name)   
        LemmeGo(keyword->value)       
        LemmeGo(keyword->pre_comment) 
        LemmeGo(keyword->line_comment)            
        LemmeGo(keyword)
    }               
                    
    return(next_kwd);
                    
}  /*  End:  "OdlFreeKwd"  */
                    
        


/*========================================================================*/
/*                                                                        */
/*                    File and File Name routines                         */
/*                                                                        */
/*========================================================================*/

/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlOpenMessageFile                                              */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a pointer to an opened file based on       */
/*      what is passed in.  If message_fptr is not NULL, the we assume  */
/*      that the file is already open and return message_fptr.  If      */
/*      message_fname is NULL, or we can't open message_fname, then     */
/*      we return stdout.  If message_fname can be opened, then we      */
/*      return a pointer to the newly opened file.                      */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

FILE *OdlOpenMessageFile (message_fname, message_fptr, suppress_messages)

const char *message_fname;
FILE *message_fptr;
int suppress_messages;

#else

FILE *OdlOpenMessageFile (const char *message_fname, FILE *message_fptr,
                          int suppress_messages)

#endif
{
    FILE *fptr = {stdout};

    if (message_fptr != NULL)
        fptr = message_fptr;
    else
        if (message_fname != NULL && ! suppress_messages)
        {
            if ((fptr = (FILE *) fopen(message_fname, "a")) == NULL)
            {
                fptr = stdout;
                OdlPrintMessage(NULL, NULL, 0, "Unable to open the output "
                    "file.  Messages will be written to the terminal",
                    suppress_messages);
            }
        }

    return(fptr);

}  /*  End routine:  "OdlOpenMessageFile"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetFileName                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine extracts the file name from a "^" keyword,         */
/*      allocates storage for it, and returns a pointer to this new     */
/*      character string.  It also returns information about where      */
/*      the data actually begins in the file.                           */
/*                                                                      */
/*      For example, lets say we're looking at a label in a file        */
/*      called test.lbl, and we want to get the file name assocated     */
/*      the FNAME keyword.  Here are the possible values this keyword   */
/*      might have, and what information would be returned for each     */
/*      possibility:                                                    */
/*                                                                      */
/*            ^FNAME = 17                                               */
/*                                                                      */
/*               file name            : test.lbl (attached)             */
/*               *start_location      : 17                              */
/*               *start_location type : ODL_RECORD_LOCATION             */
/*                                                                      */
/*          ^FNAME = 29 <RECORD>                                        */
/*                                                                      */
/*             file name            : test.lbl (attached)               */
/*             *start_location      : 29                                */
/*             *start_location type : ODL_RECORD_LOCATION               */
/*                                                                      */
/*          ^FNAME = 197 <RECORDS>                                      */
/*                                                                      */
/*             file name            : test.lbl (attached)               */
/*             *start_location      : 197                               */
/*             *start_location type : ODL_RECORD_LOCATION               */
/*                                                                      */
/*          ^FNAME = 346 <BYTE>                                         */
/*                                                                      */
/*             file name            : test.lbl (attached)               */
/*             *start_location      : 346                               */
/*             *start_location type : ODL_BYTE_LOCATION                 */
/*                                                                      */
/*          ^FNAME = 2189 <BYTES>                                       */
/*                                                                      */
/*             file name            : test.lbl (detached)               */
/*             *start_location      : 2189                              */
/*             *start_location type : ODL_BYTE_LOCATION                 */
/*                                                                      */
/*          ^FNAME = "file_name.dat"                                    */
/*                                                                      */
/*             file name            : file_name.dat (detached)          */
/*             *start_location      : 1                                 */
/*             *start_location type : ODL_RECORD_LOCATION               */
/*                                                                      */
/*          ^FNAME = ("file_name.dat", 17)                              */
/*                                                                      */
/*             file name            : file_name.dat (detached)          */
/*             *start_location      : 17                                */
/*             *start_location type : ODL_RECORD_LOCATION               */
/*                                                                      */
/*          ^FNAME = ("file_name.dat", 29 <RECORD>)                     */
/*                                                                      */
/*             file name            : file_name.dat (detached)          */
/*             *start_location      : 29                                */
/*             *start_location type : ODL_RECORD_LOCATION               */
/*                                                                      */
/*          ^FNAME = ("file_name.dat", 197 <RECORDS>)                   */
/*                                                                      */
/*             file name            : file_name.dat (detached)          */
/*             *start_location      : 197                               */
/*             *start_location type : ODL_RECORD_LOCATION               */
/*                                                                      */
/*          ^FNAME = ("file_name.dat", 346 <BYTE>)                      */
/*                                                                      */
/*             file name            : file_name.dat (detached)          */
/*             *start_location      : 346                               */
/*             *start_location type : ODL_BYTE_LOCATION                 */
/*                                                                      */
/*          ^FNAME = ("file_name.dat", 2189 <BYTES>)                    */
/*                                                                      */
/*             file name            : file_name.dat (detached)          */
/*             *start_location      : 2189                              */
/*             *start_location type : ODL_BYTE_LOCATION                 */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlGetFileName (keyword, start_location, start_location_type)

KEYWORD *keyword;
unsigned long *start_location;
unsigned short *start_location_type;

#else

char *OdlGetFileName (KEYWORD *keyword, unsigned long *start_location, 
                      unsigned short *start_location_type)

#endif
{
    char *fname = {NULL};
    char *text = {NULL};
    char *unit = {NULL};
    char *first_word = {NULL};
    char *second_word = {NULL};

    if (keyword != NULL)
    {
        /*  Make a copy of the keyword's value  */
        CopyString(text, keyword->value)
    
        /*  Get rid of parens, braces, and commas  */
        ReplaceChar(text, '(', ' ')
        ReplaceChar(text, ')', ' ')
        ReplaceChar(text, '{', ' ')
        ReplaceChar(text, '}', ' ')
        ReplaceChar(text, ',', ' ')
    
        /*  Locate the units string  */
        unit = (char *) strchr(text, '<');
    
        /*  Remove the units string if it's there  */
        if (unit != NULL) *unit = '\0';
    
        /*  Find the first word  */
        first_word = (char *) OdlFirstWord(text);
    
        /*  If the first word is quoted, then it's a file name  */
        if ((*first_word == '"') || (*first_word == '\''))
        {
            /*  Look for a second word  */
            second_word = (char *) OdlNextWord(first_word);
    
            /*  If we can't find one, then the location is record 1  */
            if (*second_word == '\0')
                *start_location = 1;
            else
            {
                /*  Otherwise, the second word is the location  */
                *start_location = atoi(second_word);
                *(second_word - 1) = '\0';
            }
    
            /*  Copy and clean up the file name  */
            CopyString(fname, (first_word+1))
            ReplaceChar(fname, '"', ' ');
            ReplaceChar(fname, '\'', ' ');
            StripTrailing(fname, ' ')
        }
        else
        {
            /*  Since the first word isn't quoted, we assume that it's a     */
            /*  location, and that the file name is the one associated with  */
            /*  the keyword itself (e.g., we're looking at attached data)    */
            *start_location = atoi(first_word);
            CopyString(fname, keyword->file_name)

        }  /*  End:  "if ((*first_word == '"') || ... else ..."  */
    
        /*  No unit string means a record location  */
        if (unit == NULL)
            *start_location_type = ODL_RECORD_LOCATION;
        else
        {
            /*  Otherwise, find out what kind of units string we have  */
            UpperCase(unit)
            *unit = '<';  /* Bug fix SM 10/24/94 */
            if (strncmp(unit, "<BYTE", 5) == 0)
                *start_location_type = ODL_BYTE_LOCATION;
            else
                *start_location_type = ODL_RECORD_LOCATION;
        }
    
        LemmeGo(text)

    }  /*  End:  "if (keyword != NULL) ..."  */

    return(fname);

}  /*  End:  "OdlGetFileName"  */




/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlGetFileSpec                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine is still TBD.  It's supposed to locate a file      */
/*      (whose name is passed in) by using the file search rules in the */
/*      PDS Standards Reference.  At the moment it just looks in the    */
/*      current directory.                                              */
/*                                                                      */
/*      WARNING:  This routine allocates memory to hold the file spec   */
/*                that is returned.  This memory will have to be freed. */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlGetFileSpec (fname)   /* this is still TBD */

char *fname;

#else

char *OdlGetFileSpec (char *fname)   /* this is still TBD */

#endif
{
    char *fspec = {NULL};

    if (fname != NULL)
    {
        /* this is still TBD */
        CopyString(fspec, fname)
    }

    return(fspec);

}  /*  End:  "OdlGetFileSpec"  */





/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlLocateStart                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine opens the filespec passed in, and attempts to      */
/*      find the start of the data by using the start_location and      */
/*      location_type that were passed in.  It returns a pointer to     */
/*      the start of the data within the file.                          */
/*                                                                      */	
/************************************************************************/

#ifdef _NO_PROTO

FILE *OdlLocateStart (filespec, start_location, start_location_type)

char *filespec;
unsigned long start_location;
unsigned short start_location_type;

#else

FILE *OdlLocateStart (char *filespec, unsigned long start_location, 
                      unsigned short start_location_type)

#endif
{
    FILE *fptr = {NULL};
    unsigned short reached_the_end = {FALSE};
    char buffer [TB_MAX_BUFFER];  /* Bug fix 11/2/94 SM:                     */
                                  /* Was TB_MAX_BUFFER + 1, which won't      */
                                  /* compile on platforms with 2-byte ints,  */
                                  /* because it's 1 bigger than than MAX_INT.*/
    unsigned long i;

    if (filespec != NULL)
    {
        if (start_location_type == ODL_BYTE_LOCATION)
        {
            fptr = (FILE *) fopen(filespec, "rb");
            if ((fptr != NULL) && (start_location > 1))
                reached_the_end = (fseek(fptr,start_location,0) != 0);
        }
        else
        {
            fptr = (FILE *) fopen(filespec, "r");
            if (fptr != NULL)
            {
                for (i=1; ((i < start_location) && (! reached_the_end)); ++i)
                {
                    if (! fgets(buffer, TB_MAX_BUFFER, fptr))
                        reached_the_end = TRUE;
                }
            }

        }  /*  End:  "if (start_location_type == ... else ..."  */
 
        if (reached_the_end) CloseMe(fptr)

    }  /*  End:  "if (filespec != NULL) ..."  */

    return(fptr);

}  /*  End:  "OdlGetFileSpec"  */






/*========================================================================*/
/*                                                                        */
/*                            Print Routines                              */
/*                                                                        */
/*========================================================================*/

/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPrintMessage                                                 */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine prints a formatted message either to stdout or to  */
/*      a message file, depending on what was passed in.  Messages are  */
/*      formatted to look like this:                                    */
/*                                                                      */
/*         <line number> -- <message text>                              */
/*                                                                      */
/*      If the line_number is zero, then just the message text is       */
/*      printed, with no formatting.                                    */
/*                                                                      */
/************************************************************************/
char ODLErrorMessage[120] = "";

#ifdef _NO_PROTO

short OdlPrintMessage (message_fname, message_fptr, line_number, text,
                       suppress_messages)

const char *message_fname;
FILE *message_fptr;
long line_number;
const char *text;
int suppress_messages;

#else

short OdlPrintMessage (const char *message_fname, FILE *message_fptr, 
                       long line_number, const char *text, 
                       int suppress_messages)

#endif
{
    FILE *m_ptr = {NULL};
    char line_prompt[20];
    char *line_out = {NULL};

    ++odl_message_count;

    if (!suppress_messages)
    {
        m_ptr = (message_fptr != NULL) ? message_fptr :
                 (FILE *) OdlOpenMessageFile(message_fname, message_fptr,
                                             suppress_messages);
    
        if (line_number == 0)
            (void)strcpy(line_prompt, "");
        else
            (void)sprintf(line_prompt, " Line %ld -- ", line_number);
    
        if (text == NULL)
        {
            NewString(line_out, (20 + (long)strlen(line_prompt)))
            (void)sprintf(line_out, "%s Unknown error", line_prompt);
        }
        else
        {
            NewString(line_out, (20 + (long)strlen(line_prompt) + (long)strlen(text)))
            (void)sprintf(line_out, "%s%s", line_prompt, text);
        }
    
        line_out = OdlFormatMessage(line_out);
        (void)strncpy(ODLErrorMessage, line_out, sizeof(ODLErrorMessage));
        ODLErrorMessage[sizeof(ODLErrorMessage) - 1] = '\0';
        (void)fprintf(m_ptr, "%s", line_out);
        LemmeGo(line_out)
    
        /*  if we opened the message file in this routine then close it  */
        if ((m_ptr != stdout) && (message_fptr == NULL))
            CloseMe(m_ptr)

    }  /*  End:  "if (!suppress_messages) ..."  */
    else
    {
        /* Copy the message to the error string, but only if a message isn't
           already there */
        if (strcmp(ODLErrorMessage, "") == 0)
        {
            (void)strncpy(ODLErrorMessage, text, sizeof(ODLErrorMessage));
            ODLErrorMessage[sizeof(ODLErrorMessage) - 1] = '\0';
        }
    }
    
    return(FALSE);

}  /*  End routine:  OdlPrintMessage  */

/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPrintLine                                                    */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      J. S. Hughes (Jet Propulsion Laboratory)                        */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    February 20, 1995                                        */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*	02-20-95    jsh - a copy of OdlPrintMessage cleaned up for      */
/*		    simple output of lines to message file		*/
/*                  NOTE: '\n' is not appended or assumed               */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine prints a simple  line either to stdout or to 	*/
/*      a message file, depending on what was passed in.                */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

short OdlPrintLine (message_fname, message_fptr, text, suppress_messages)

const char *message_fname;
FILE *message_fptr;
const char *text;
int suppress_messages;

#else

short OdlPrintLine (const char *message_fname, FILE *message_fptr, 
                       const char *text, int suppress_messages)

#endif
{
    FILE *m_ptr = {NULL};

    if (!suppress_messages)
    {
        m_ptr = (message_fptr != NULL) ? message_fptr :
                 (FILE *) OdlOpenMessageFile(message_fname, message_fptr,
                                             suppress_messages);
    
        if (text == NULL)
        {
	    (void)fprintf(m_ptr, "%s", "Unknown error\n");
        }
        else
        {
            (void)fprintf(m_ptr, "%s", text);
        }
    
        /*  if we opened the message file in this routine then close it  */
        if ((m_ptr != stdout) && (message_fptr == NULL))
            CloseMe(m_ptr)

    }  /*  End:  "if (!suppress_messages) ..."  */
    
    return(FALSE);

}  /*  End routine:  OdlPrintLine  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

char *OdlFormatMessage (text)
                                                      
char *text;

#else

char *OdlFormatMessage (char *text)

#endif
{
    char *new_text = {NULL};
    char *first_char = {NULL};
    char *last_char = {NULL};
    char *dashes = {NULL};
    char *blanks = {NULL};
    char *c = {NULL};
    char save_it = {'\0'};
    long report_indent = {0};
    long report_width = {75};
    long line_size = {0};
    long len = {0};
    long i = {0};

    /* IF a message was passed in THEN                                     */
    if (text != NULL)
    {
        NewString(new_text, 1)

        /* Find the double dash delimiter thingy in the message.           */
        /*     Messages will look something like this:                     */
        /*         WARNING: Line 123 -- BANDS: Not in data dictionary.     */
        /*     We are using the location of the " -- " characters to       */
        /*     figure out how far the wrapped part of the line should      */
        /*     be indented.                                                */

        if ((dashes = strstr(text, " -- ")) != NULL)
            report_indent = 4 + ((long) (dashes - text));

        if (report_indent >= (report_width - 2))
            report_indent = 0;

        /* Initialize the string of blanks used for indentation.           */
	NewString(blanks, report_indent + 1)
        for (i=0; i < report_indent; ++i)
	    *(blanks+i) = ' ';
	*(blanks+i) = '\0';

        /* Figure out the size of the wrapped parts of the line.           */
        line_size = report_width - report_indent;

        /* Now that we have all that out of the way, we can LOOP through   */
        /*         the string until we have wrapped and written the        */
        /*         whole thing.                                            */
        for (first_char=text; *first_char != '\0'; first_char = last_char)
        {
            /* Find the length of the remaining part of the string.        */
            len = strlen(first_char);

            /* IF we are at the beginning of the string THEN               */
            /*     Use the total width of the report to figure out where   */
            /*         the end of the line should be.                      */
            /* ELSE                                                        */
            /*     Write the blanks to the report file and use the space   */
            /*         left over after indentation to figure out where     */
            /*         the end of the line should be.                      */
            /* ENDIF                                                       */

            if (first_char == text) 
            {
                if (len > report_width)
		    last_char = (char *) (first_char + report_width);
                else
		    last_char = (char *) (first_char + len);
            }
            else
            {
                AppendString(new_text, blanks)

                if (len > line_size)
		    last_char = (char *) (first_char + line_size);
                else
		    last_char = (char *) (first_char + len);

            }  /*  End:  "if (first_char == text) ... else ..."  */

            /* IF the current part of the message is still too large to    */
            /*        fit without wrapping THEN                            */
            /*     Find the last blank in the line and wrap there.         */
            /* ENDIF                                                       */

            if (*last_char != '\0')
            {
                for (c = last_char; ((c >= first_char) && (*c != ' ')); --c) ;
                           
                if (c > first_char)
                    last_char = c;

            }  /*  End: "if (*last_char != '\0') ..."  */

            /* Append the current part of the message onto the new string  */
            save_it = *last_char;
            *last_char = '\0';
            AppendString(new_text, first_char)
            AppendString(new_text, "\n")
            *last_char = save_it;

            /* Bypass the last blank character.                            */
            if (*last_char == ' ')
                ++last_char;

        }  /*  End:  "for (first_char = text; ..."  */

        /* Deallocate local storage.                                       */
        LemmeGo(blanks)

    }  /*  End:  "if ((text != NULL) && ..."  */
       
    LemmeGo(text)

    return(new_text);

}  /*  End routine:  "OdlFormatMessage"  */




/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPrintHierarchy                                               */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine prints the object hierarchy to a message file.     */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

void OdlPrintHierarchy (object, message_fname, message_fptr, suppress_messages)

OBJDESC *object;
char *message_fname;
FILE *message_fptr;
int suppress_messages;

#else

void OdlPrintHierarchy (OBJDESC *object, char *message_fname, 
                        FILE *message_fptr, int suppress_messages)

#endif
{
    OBJDESC *obj = {NULL};
    KEYWORD *kwd = {NULL};
    FILE *m_ptr = {NULL};
    char *format = {NULL};
    const char *no_name = {"<no name>"};
    char msgtext [TB_MAXLINE + 1];

    m_ptr = (message_fptr != NULL) ? message_fptr :
                 (FILE *) OdlOpenMessageFile(message_fname, message_fptr,
                                             suppress_messages);

    for (obj=object; obj != NULL; obj=(OBJDESC *) OdlTraverseTree(obj, object->level))
    {
        NewString(format, (TB_MAXLINE + 1))

        kwd = (KEYWORD *)OdlFindKwd(obj, "NAME", (char *)NULL,(long)1, ODL_THIS_OBJECT);

        if ((kwd == NULL) || (kwd->value == NULL))
        {
            (void)sprintf(format, " Line %-5lu %%%lud %%%lus", 
                    obj->line_number, (obj->level + 1),
                    (2*(obj->level) + strlen(no_name)));

            (void)sprintf(msgtext, format, obj->level, no_name);
            OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        }
        else
        {
            (void)sprintf(format, " Line %-5lu %%%lud %%%lus", 
                    obj->line_number, (obj->level + 1),
                    (2*(obj->level) + strlen(kwd->value)));
            (void)sprintf(msgtext, format, obj->level, kwd->value);
            OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        }

        if (obj->class == NULL)
        {
            OdlPrintLine(message_fname, m_ptr, "  --  <no class>\n",
                         suppress_messages);
        }
        else {
            (void)sprintf(msgtext, "  --  %s\n", obj->class);
            OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
	}

        LemmeGo(format)

    }  /*  End:  "for (obj=object; ..."  */

    /*  if we opened the message file in this routine then close it  */
    if ((m_ptr != stdout) && (message_fptr == NULL))
        CloseMe(m_ptr)
    
    return;

}  /*  End routine:  "OdlPrintHierarchy"  */




/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlPrintLabel                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*      03-11-97    Add GROUP logic                                     */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine prints the ODL tree to a message file, in ODL      */
/*      format, unless the suppress_messages flag is set.               */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

void OdlPrintLabel (object, message_fname, message_fptr, root_level)

OBJDESC *object;
char *message_fname;
FILE *message_fptr;
unsigned long root_level;

#else

void OdlPrintLabel (OBJDESC *object, char *message_fname, FILE *message_fptr, 
                    unsigned long root_level)

#endif
{
    FILE *m_ptr = {NULL};
    char *blanks = {NULL};
    int i;
    char msgtext [TB_MAXLINE + 1];
    /* declare suppress messages as a local variable in case it never needs to
       be changed into a parameter */
    int suppress_messages = 0;

    if (!suppress_messages)
    {
        m_ptr = (message_fptr != NULL) ? message_fptr :
                     (FILE *) OdlOpenMessageFile(message_fname, message_fptr,
                                                 suppress_messages);
        if (object != NULL)
        {
            NewString(blanks, (4*object->level))
            for (i=1; i < object->level; ++i) (void)strcat(blanks, "  ");
    
            if (object->pre_comment != NULL)
            {
                OdlPrintLine(message_fname, m_ptr, object->pre_comment,
                                   suppress_messages);
            }
    
            if (object->parent != NULL)
            {
	        if (object->is_a_group == ODL_OBJECT)
		{
		    if (object->class == NULL) {
                (void)sprintf(msgtext, "%sOBJECT", blanks);
                OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
		    }
            else {
		        (void)sprintf(msgtext, "%sOBJECT = %s", blanks, object->class);
		        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
	   	    }
		}
	        else
		{
		    if (object->class == NULL) {
                        (void)sprintf(msgtext, "%sGROUP", blanks);
                        OdlPrintLine(message_fname, m_ptr, msgtext,
                                     suppress_messages);
	  	    }
                    else {
		        (void)sprintf(msgtext, "%sGROUP = %s", blanks, object->class);
		        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
	  	    }
		}

	        if ((object->line_comment != NULL) 
                && (*object->line_comment != '\0')) {
                OdlPrintLine(message_fname, m_ptr, " ", suppress_messages);
                OdlPrintLine(message_fname, m_ptr, object->line_comment, 
                             suppress_messages);
		}
    
                OdlPrintLine(message_fname, m_ptr, "\n", suppress_messages);
    
            }  /*  End:  "if (object->parent != NULL) ..."  */
    
            OdlPrintKeywords(object, NULL, m_ptr, suppress_messages);
            OdlPrintLabel(object->first_child, (char *)NULL, m_ptr, root_level);
    
            if (object->post_comment != NULL)
            {
                OdlPrintLine(message_fname, m_ptr, object->post_comment,
                             suppress_messages);
            }
    
            if (object->parent != NULL)
            {
	        if (object->is_a_group == ODL_OBJECT)
		{
                    if (object->class == NULL) {
                        (void)sprintf(msgtext, "%sEND_OBJECT", blanks);
                        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
	   	    }
                    else {
                        (void)sprintf(msgtext, "%sEND_OBJECT = %s", blanks, object->class);
                        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
		    }
  		}
		else
		{
                    if (object->class == NULL) {
                        (void)sprintf(msgtext, "%sEND_GROUP", blanks);
                        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
	 	    }
                    else {
                        (void)sprintf(msgtext, "%sEND_GROUP = %s", blanks, object->class);
                        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
		    }
      		}

                if ((object->end_comment != NULL) 
                    && (*object->end_comment != '\0')) {
                    OdlPrintLine(message_fname, m_ptr, " ", suppress_messages);
                    OdlPrintLine(message_fname, m_ptr, object->end_comment, suppress_messages);
		}

                OdlPrintLine(message_fname, m_ptr, "\n", suppress_messages);
    
            }  /*  End:  "if (object->parent != NULL) ..."  */
    
            if (object->level > root_level)
                OdlPrintLabel(object->right_sibling, (char *)0, m_ptr, root_level);
    
            LemmeGo(blanks)
    
            if (object->parent == NULL)
                OdlPrintLine(message_fname, m_ptr, "END\n", suppress_messages);
    
        }  /*  End:  "if (object != NULL) ..."  */
    
        /*  if we opened the message file in this routine then close it  */
        if ((m_ptr != stdout) && (message_fptr == NULL))
            CloseMe(m_ptr)
    
    }  /*  End:  "if (!suppress_messages) ..."  */

    return;

}  /*  End routine:  "OdlPrintLabel"  */



/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

void OdlPrintKeywords (object, message_fname, message_fptr, suppress_messages)

OBJDESC *object;
char *message_fname;
FILE *message_fptr;
int suppress_messages;

#else

void OdlPrintKeywords (OBJDESC *object, char *message_fname, 
                              FILE *message_fptr, int suppress_messages)

#endif
{
    KEYWORD *keyword = {NULL};
    FILE *m_ptr = {NULL};
    short sfdu_only = {FALSE};
    char *blanks = {NULL};
    int i;

    m_ptr = (message_fptr != NULL) ? message_fptr :
                 (FILE *) OdlOpenMessageFile(message_fname, message_fptr,
                                             suppress_messages);

    if (object != NULL)
    {
        NewString(blanks, (4*object->level))
        if (object->level > 0)
            for (i=0; i < object->level; ++i) (void)strcat(blanks, "  ");

        for (keyword=object->first_keyword; keyword != NULL; 
                 keyword = keyword->right_sibling)
        {
            if (keyword->pre_comment != NULL)
                OdlPrintLine(message_fname, m_ptr, keyword->pre_comment, suppress_messages);
 
            OdlPrintLine(message_fname, m_ptr, blanks, suppress_messages);

            sfdu_only = FALSE;
            if (keyword->name == NULL)
                OdlPrintLine(message_fname, m_ptr, "unknown_keyword", suppress_messages);
            else
            { 
                OdlPrintLine(message_fname, m_ptr, keyword->name, suppress_messages);
                sfdu_only = ((strncmp(keyword->name, "NJPL", 4) == 0) ||
                                (strncmp(keyword->name, "CCSD", 4) == 0));
 	    }

            if ((keyword->value != NULL) && (! sfdu_only)) {
                OdlPrintLine(message_fname, m_ptr, " = ", suppress_messages);
                OdlPrintLine(message_fname, m_ptr, keyword->value, suppress_messages);
            }

            if ((keyword->line_comment != NULL) 
                && (*keyword->line_comment != '\0')) {
                OdlPrintLine(message_fname, m_ptr, " ", suppress_messages);
                OdlPrintLine(message_fname, m_ptr, keyword->line_comment, suppress_messages);
            }

            OdlPrintLine(message_fname, m_ptr, "\n", suppress_messages);

        }  /*  End:  "for (keyword=object ..."  */

        LemmeGo(blanks)

    }  /*  End:  "if (object != NULL) ..."  */

    /*  if we opened the message file in this routine then close it  */
    if ((m_ptr != stdout) && (message_fptr == NULL))
        CloseMe(m_ptr)
    
    return;

}  /*  End routine:  "OdlPrintKeywords"  */




/*========================================================================*/
/*                                                                        */
/*                       Parser-specific routines                         */
/*                                                                        */
/*========================================================================*/


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlParseFile                                                    */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*      03-11/97    Add GROUPS objects                                  */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine actually does the parsing of a label file and      */
/*      returns a pointer to the root object of the tree.               */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlParseFile(label_fname, label_fptr, message_fname, message_fptr,
                     suppress_messages, suppress_metrics, suppress_hierarchy,
                     ignore_missing_end)

char *label_fname;
FILE *label_fptr;
char *message_fname;
FILE *message_fptr;
int suppress_messages;
unsigned short suppress_metrics;
unsigned short suppress_hierarchy;
unsigned short ignore_missing_end;

#else

OBJDESC *OdlParseFile( char *label_fname, FILE *label_fptr, 
                      char *message_fname, FILE *message_fptr,
                      int suppress_messages, 
                      unsigned short suppress_metrics, 
                      unsigned short suppress_hierarchy,
                      unsigned short ignore_missing_end)

#endif        

{
    OBJDESC *root = {NULL};
    OBJDESC *curr_object = {NULL};
    KEYWORD *curr_keyword = {NULL};
    FILE *m_ptr = {NULL};
    FILE *l_ptr = {NULL};
    char *left_part = {NULL};
    char *equals = {NULL};
    char *right_part = {NULL};
    char *comment = {NULL};
    char *c = {NULL};
    char *tc = {NULL};
    char *tintext = {NULL};
    char *text = {NULL};
    char line_comment [TB_MAXLINE + 1];
    char intext [TB_MAXLINE + 1];
    char msgtext [TB_MAXLINE + 1];
    long line_number = {0};
    long value_list_line_number = {0};
    long object_count = {0};
    long end_object_count = {0};
    long keyword_count = {0};
    long comment_count = {0};
    long brace_nesting = {0};
    long paren_nesting = {0};
    short end_found = {FALSE};
    short value_list = {FALSE};
    short equals_found = {FALSE};
    short val_found = {FALSE};
    short balanced = {FALSE};
    short oddquotes = {FALSE};
    unsigned short is_a_group;

    odl_message_count = 0;

    /*  either use the file pointer passed in or open the message file  */
    m_ptr = (message_fptr != NULL) ? message_fptr :
                 (FILE *) OdlOpenMessageFile(message_fname, message_fptr,
                                             suppress_messages);
    /*  opening remarks  */
    if (label_fname == NULL)
        (void)sprintf(msgtext, "Parsing File:  (no file name provided)");
    else
        (void)sprintf(msgtext, "Parsing File:  %s", label_fname);

    OdlPrintLine(message_fname, m_ptr, "\n--------------------------------------------------------------------------\n", suppress_messages);
    OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
    OdlPrintLine(message_fname, m_ptr, "\n--------------------------------------------------------------------------\n\n", suppress_messages);

    /*  either use the file pointer passed in or open the label file  */
    l_ptr = (label_fptr != NULL) ? label_fptr :
                                   (label_fname == NULL) ? NULL :
                                                           (FILE *) fopen(label_fname,"r");
    if (l_ptr == NULL)
    {
        (void)OdlPrintMessage(message_fname, m_ptr, (long)0,
            "Unable to open the label file.  Parsing cannot continue", suppress_messages);
    }
    else
    {
        NewString(comment, 1)

        /*  Initialize a ROOT object to start the tree  */
        curr_object = root = OdlNewObjDesc("ROOT",(char *)0,(char *)0,(char *)0,(char *)0,label_fname,0,(long)0);

        /*  read the label file  */
        while (! end_found && fgets(intext, TB_MAXLINE, l_ptr))
        {
            ++line_number;

            StripUnprintables(intext)  /*  removes linefeeds and such  */
            ReplaceChar(intext, '	', ' ')  /*  turns TABs into blanks  */
            StripTrailing(intext, ' ') /*  removes trailing blanks     */

            /* locate and skip SFDU */

            if (line_number == 1L)
            {
               if (! strncmp(intext, "CCSD", 4))
                   continue;
            }

            /*  locate, save, and remove comment text from the line  */
            *line_comment = '\0';
            tintext = intext;
            oddquotes = FALSE;
            while ((c = strstr(tintext, "/*")) != NULL)
            {
		for (tc = tintext; tc < c; tc++)
			if (*tc == '"')
            {
				if (oddquotes) oddquotes = FALSE;
				else oddquotes = TRUE;
            }
		if (! oddquotes)  
                {
			oddquotes = FALSE;
			++comment_count;
			(void)strcpy(line_comment, c);
			*c = '\0';
			StripTrailing(tintext, ' ');
			break; 
		}
		else
			tintext = c+1;
            }

            c = OdlFirstWord(intext);

            if (text != NULL || *c != '\0')
            {
                AppendString(text, intext)                

                c = OdlFirstWord(text);

                if (strcmp(c, "END") == 0)
                {
                    balanced = TRUE;
                    end_found = TRUE;
                    break;
                }
                else
                if (strcmp(c, "END_OBJECT") == 0 || strcmp(c, "END_GROUP") == 0 )
                {
                    balanced = TRUE;
                }
                else
                {
                    if (! equals_found)
                        equals_found = (strchr(text, '=') != NULL);
    
                    if (! val_found && equals_found)
                    {
                        c = (char *) LastChar(text);
                        val_found = (*c != '=');
                    }

                    if (val_found && equals_found)
                        balanced = CheckBalance(text);

                    if (! balanced)
  	                AppendString(text, "\n")                
                }
	    }

            if (balanced)
            {
                /*  locate the keyword, the equals sign, and the value  */
                left_part = OdlFirstWord(text);

                if ((equals = (char *)strchr(left_part, '=')) != NULL)
                    right_part = OdlFirstWord(equals+1);
                else
                {
                    equals = text + strlen(text);
                    right_part = equals;
                }

        /*------------------------------------------------------------------*/
        /*  Here's where the parsing begins.  First, we take care of three  */
        /*  special cases:  multi-line quoted values, multi-line value      */
        /*  lists, and blank lines.  If the current line isn`t one          */
        /*  of these than it's either an OBJECT statement, an END_OBJECT    */
        /*  statement, the END of the label, or a new KEYWORD.              */
        /*------------------------------------------------------------------*/

                /*  we've discovered the beginning of a new object  */
                if ((strncmp(left_part, "OBJECT ", 7) == 0) || (strcmp(left_part, "OBJECT") == 0) ||
                    (strncmp(left_part, "GROUP ", 6) == 0) || (strcmp(left_part, "GROUP") == 0))

                {
                    ++object_count;
                    ++(curr_object->child_count);

		    if ((strncmp(left_part, "OBJECT ", 7) == 0) || (strcmp(left_part, "OBJECT") == 0))
		      is_a_group = ODL_OBJECT;
		    else
		      is_a_group = ODL_GROUP;
    
                    /*  validate the new object's class identifier  */
                    (void)OdlValidObjDesc(curr_object, equals, right_part,
                                     message_fname, m_ptr, line_number, suppress_messages);
    
                    /*  make the new object a child of the current object  */
                    curr_object = OdlPasteObjDesc(OdlNewObjDesc(right_part,
                                                  comment,line_comment,(char *)0,(char *)0,
                                                  label_fname,is_a_group,line_number),
                                                curr_object);
    
                    /*  reset the comment text string  */
                    LemmeGo(comment)
                    NewString(comment, 1)
                }
                else
        /*------------------------------------------------------------------*/
                /*  we've discovered the end of the current object  */
                if ((strncmp(left_part, "END_OBJECT ", 11) == 0) || (strcmp(left_part, "END_OBJECT") == 0) ||
                    (strncmp(left_part, "END_GROUP ", 10) == 0) || (strcmp(left_part, "END_GROUP") == 0))
                {
                    ++end_object_count;

                    if ((strncmp(left_part, "END_OBJECT ", 11) == 0) || (strcmp(left_part, "END_OBJECT") == 0))
		      is_a_group = ODL_OBJECT;
		    else
		      is_a_group = ODL_GROUP;
    
                    /*  validate the end_object's class identifier  */
                    (void)OdlValidEndObjDesc(curr_object, equals, right_part,
                                   message_fname, m_ptr, line_number, is_a_group, 
                                   suppress_messages);
    
                    /*  set the current object's remaining comment fields  */
                    CopyString(curr_object->post_comment, comment)
                    CopyString(curr_object->end_comment, line_comment)
    
                    /*  make curr object's parent the new current object  */
                    if (curr_object->parent != NULL)
                        curr_object = curr_object->parent;
    
                    /*  reset the comment text string  */
                    LemmeGo(comment)
                    NewString(comment, 1)
                }
                else
        /*------------------------------------------------------------------*/
                /*  we've reached the end of the label  */
                if ((strncmp(left_part, "END ", 4) == 0) || 
                     (strcmp(left_part, "END") == 0))
                {
                    end_found = TRUE;
                    CopyString(curr_object->post_comment, comment)
                }
                else
        /*------------------------------------------------------------------*/
                /*  We've discovered a keyword and its value  */
                {
                    ++keyword_count;
    
                    /*  validate the keyword and its values  */
                    (void)OdlValidKwd(curr_object, left_part, equals, 
                             right_part, message_fname, m_ptr, line_number, suppress_messages);
    
                    /*  Add the keyword to the current object  */
                    curr_keyword = OdlPasteKwd(OdlNewKwd(left_part, right_part,
                                                     comment, line_comment, 
                                                     label_fname, line_number),
                                                  curr_object);
    
                    /*  we've got a potential multi-line value list if the 
                        first character of the value is either an open 
                        brace, '{', or an open paren, '('
                    */

                    if ((value_list = curr_keyword->is_a_list) == TRUE)
                    {
                        /*  validate that the braces and parens are correct  */
                        (void)OdlValidBraces(curr_keyword->value,
                                    brace_nesting, paren_nesting,
                                        message_fname, m_ptr, 
                                        value_list_line_number, suppress_messages);


                    /*  reset the comment text string  */
                    LemmeGo(comment)
                    NewString(comment, 1)
    
                }  /*  End:  "if ((strncmp(left_part, ... else ... else ..." */
    /*------------------------------------------------------------------*/
                }  /*  End:  "if (quoted_value) ... else ... else ..."  */

                equals_found = FALSE;
                val_found = FALSE;
                balanced = FALSE;
                LemmeGo(text)

	    }
    /*------------------------------------------------------------------*/
        }  /*  End:  "while (fgets(text, ..."  */

        /* if we're not sitting at the root then not enough END_OBJECTs found */
        if (curr_object->parent != NULL)
        {
            (void)OdlPrintMessage(message_fname, m_ptr, line_number,
                "Not enough END_OBJECT statements.  Some objects may be incomplete",
                suppress_messages);
        }

        /*  hey, we didn't find an end statement!  */
        if ((! end_found) && (! ignore_missing_end))
        {
            (void)OdlPrintMessage(message_fname, m_ptr, line_number,
                "END statement is missing", suppress_messages);
        }

        /*  oops, there was nothing in the label file to parse  */
        if (line_number == 0) 
            root = OdlFreeTree(root);

        LemmeGo(comment)

    }  /*  End:  "if (l_ptr == NULL) ... else ..."  */

    /*  how'd we do?  */
    if (! suppress_metrics && ! suppress_messages) 
    {
        OdlPrintLine(message_fname, m_ptr, "\n", suppress_messages);
        OdlPrintLine(message_fname, m_ptr, "           |-------------------------------|\n", suppress_messages);
        OdlPrintLine(message_fname, m_ptr, "           | Parsing Metrics:              |\n", suppress_messages);
        OdlPrintLine(message_fname, m_ptr, "           |                               |\n", suppress_messages);
        (void)sprintf(msgtext, "           | %7ld Syntax Messages       |\n", odl_message_count);
        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        OdlPrintLine(message_fname, m_ptr, "           |                               |\n", suppress_messages);
        (void)sprintf(msgtext, "           | %7ld OBJECT Statements     |\n", object_count);
        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        (void)sprintf(msgtext, "           | %7ld END_OBJECT Statements |\n", end_object_count);
        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        (void)sprintf(msgtext, "           | %7ld Keywords              |\n", keyword_count);
        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        (void)sprintf(msgtext, "           | %7ld Comments              |\n", comment_count);
        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        OdlPrintLine(message_fname, m_ptr, "           |-------------------------------|\n\n", suppress_messages);

    }  /*  End:  "if (! suppress_metrics) ..."  */

    /*  display the object hierarchy  */
    if (! suppress_hierarchy && ! suppress_messages)
    {
        if (label_fname == NULL)
            (void)sprintf(msgtext, "Object Hierarchy in File:  (no file name provided)");
        else
            (void)sprintf(msgtext, "Object Hierarchy in File:  %s", label_fname);

        OdlPrintLine(message_fname, m_ptr, "\n--------------------------------------------------------------------------\n", suppress_messages);
        OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
        OdlPrintLine(message_fname, m_ptr, "\n--------------------------------------------------------------------------\n\n", suppress_messages);

        OdlPrintHierarchy(root, message_fname, m_ptr, suppress_messages);

    }  /*  End:  "if (suppress_hierarchy) ..."  */

    /*  closing remarks  */



    if (label_fname == NULL)
        (void)sprintf(msgtext, "End of Parsing File:  (no file name provided)");
    else
        (void)sprintf(msgtext, "End of Parsing File:  %s", label_fname);

    OdlPrintLine(message_fname, m_ptr, "\n--------------------------------------------------------------------------\n", suppress_messages);
    OdlPrintLine(message_fname, m_ptr, msgtext, suppress_messages);
    OdlPrintLine(message_fname, m_ptr, "\n--------------------------------------------------------------------------\n\n", suppress_messages);

    /*  if we opened the label file in this routine then close it  */
    if (label_fptr == NULL)
        CloseMe(l_ptr)

    /*  if we opened the message file in this routine then close it  */
    if ((m_ptr != stdout) && (message_fptr == NULL))
        CloseMe(m_ptr)
    
    LemmeGo(text)

    return (root);

}  /*  End routine:  OdlParseFile  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlNestingLevel (text, brace_nesting, paren_nesting)
char *text;
long *brace_nesting;
long *paren_nesting;

#else

short OdlNestingLevel (char *text, long *brace_nesting, 
                              long *paren_nesting)

#endif
{
    char *c = {NULL};

    for (c=text; *c != '\0'; ++c)
    {
        if (*c == '{') 
            ++(*brace_nesting);
        else
        if (*c == '}') 
            --(*brace_nesting);
        else
        if (*c == '(') 
            ++(*paren_nesting);
        else
        if (*c == ')') 
            --(*paren_nesting);
    }

    return((*brace_nesting == 0) && (*paren_nesting == 0));

}  /*  End routine:  "OdlNestingLevel"  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlValidBraces (text, brace_nesting, paren_nesting,
                          message_fname, message_fptr, line_number,
                          suppress_messages)

char *text;
long brace_nesting;
long paren_nesting;
char *message_fname;
FILE *message_fptr;
long line_number;
int suppress_messages;

#else

short OdlValidBraces (char *text, long brace_nesting, 
                             long paren_nesting, char *message_fname, 
                             FILE *message_fptr, long line_number,
                             int suppress_messages)

#endif
{
    char *c = {NULL};
    char *sp = {NULL};
    char *nesting_stack = {NULL};
    short status = {TRUE};

    /*  allocate storage for the nesting stack  */
    NewString(nesting_stack, (long)strlen(text))

    /*  validate that all braces and parens are correctly nested  */
    for (c=text,sp=(nesting_stack-1); ((*c != '\0') && (status == TRUE)); ++c)
    {
        /*  push brace or paren onto the nesting stack  */
        if ((*c == '{') || (*c == '('))
            *(++sp) = *c;
        else
        /*  nesting is ok so far, pop the nesting stack  */
        if (((*c == '}') && (*sp == '{')) || ((*c == ')') && (*sp == '(')))
            --sp;
        else
        /*  found a right brace that doesn't have a matching left one  */
        if ((*c == '}') && (*sp != '{'))
        {
            status = OdlPrintMessage(message_fname,message_fptr,line_number, 
             "Bad nesting in VALUE LIST.  Expected a right parenthesis and found a brace instead.",
             suppress_messages);
        }
        else
        /*  found a right paren that doesn't have a matching left one  */
        if ((*c == ')') && (*sp != '('))
        {
            status = OdlPrintMessage(message_fname,message_fptr,line_number, 
             "Bad nesting in VALUE LIST.  Expected a right brace and found a parenthesis instead.",
             suppress_messages);
        }

        /*  we've reached nesting level zero before reaching the end  */
        if ((sp < nesting_stack) && (*(c+1) != '\0'))
        {
            status = OdlPrintMessage(message_fname,message_fptr,line_number, 
                         "VALUE LIST not properly enclosed in braces or parentheses",
                         suppress_messages);
        }

    }  /*  End:  "for (c=text,sp=(nesting_stack-1); ..."  */

    LemmeGo(nesting_stack)

    if (brace_nesting < 0)
    {
        status = OdlPrintMessage(message_fname, message_fptr, line_number, 
                     "Too many right braces in VALUE LIST", suppress_messages);
    }
    else
    if (brace_nesting > 0)
    {
        status = OdlPrintMessage(message_fname, message_fptr, line_number, 
                     "Too many left braces in VALUE LIST", suppress_messages);
    }

    if (paren_nesting < 0)
    {
        status = OdlPrintMessage(message_fname, message_fptr, line_number, 
                     "Too many right parentheses in VALUE LIST", suppress_messages);
    }
    else
    if (paren_nesting > 0)
    {
        status = OdlPrintMessage(message_fname, message_fptr, line_number, 
                     "Too many left parentheses in VALUE LIST", suppress_messages);
    }

    return(status);

}  /*  End routine:  "OdlValidBraces"  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlValidElement (text, message_fname, message_fptr, line_number, 
                           element_number, suppress_messages)
char *text;
char *message_fname;
FILE *message_fptr;
long line_number;
long element_number;
int suppress_messages;

#else

short OdlValidElement (char *text, char *message_fname, 
                              FILE *message_fptr, long line_number, 
                              long element_number, int suppress_messages)

#endif
{
    char *message = NULL;
    char element_prompt[TB_MAXLINE + 1];
    char *save_units = 0;
    char *first_blank = {NULL};
    char *first_char = {NULL};
    char *last_char = {NULL};
    char *units_start = {NULL};
    char *units_end = {NULL};
    char *single_quote = {NULL};
    char *double_quote = {NULL};
    short status = {TRUE};

    if (element_number <= 0)
       (void)strcpy(element_prompt, "");
    else
       (void)sprintf(element_prompt, " LIST element %ld", element_number);

    single_quote = (char *) strchr(text+1, (int) '\'');
    double_quote = (char *) strchr(text+1, (int) '"');
    first_blank  = (char *) strchr(text+1, (int) ' ');
    first_char = text;
    last_char = (char *) LastChar(text);

    NewString(message, (TB_MAXLINE+(long)strlen(text)))

    /*  double quote found in the middle of the value  */
    if ((double_quote > first_char) && (double_quote < last_char))
    {
        (void)sprintf(message, "Embedded double quote in VALUE%s", element_prompt);
        status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
    }
    else
    /*  value is double quoted - everything is okay  */
    if (*first_char == '"')
    {
        status = TRUE;
    }
    else
    /*  single quote found in the middle of the value  */
    if ((single_quote > first_char) && (single_quote < last_char))
    {
        (void)sprintf(message, "Embedded single quote in VALUE%s", element_prompt);
        status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
    }
    else
    /*  value is single quoted - fine if not just a quote  */
    if ((*first_char == '\'') && (*last_char == '\''))
    {
        if (first_char == last_char)
        {
            (void)sprintf(message, "Unpaired single quote in VALUE%s", element_prompt);
            status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
        }
    }
    else
    /*  value is missing a closing single quote  */
    if ((*first_char == '\'') && (*last_char != '\''))
    {
        (void)sprintf(message, "Unpaired single quote in VALUE%s", element_prompt);
        status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
    }
    else
    /*  value is missing an opening single quote  */
    if ((*first_char != '\'') && (*last_char == '\''))
    {
        (void)sprintf(message, "Unpaired single quote in VALUE%s", element_prompt);
        status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
    }
    else
    /*  value is missing an opening double quote  */
    if ((*first_char != '"') && (*last_char == '"'))
    {
        (void)sprintf(message, "Unpaired double quote in VALUE%s", element_prompt);
        status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
    }
    else
    /*  current value list element is just a double quote  */
    if ((element_number > 0) && 
        (first_char == (last_char-1)) && (*first_char == '"'))
    {
        (void)sprintf(message, "Unpaired double quote in VALUE%s", element_prompt);
        status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
    }
    else
    /*  current value list element is missing a closing double quote  */
    if ((element_number > 0) && 
        (*first_char == '"') && (*(last_char-1) != '"'))
    {
        (void)sprintf(message, "Unpaired double quote in VALUE%s", element_prompt);
        status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
    }
    else
    /*  value is unquoted  */
    if ((*first_char != '\'') && (*last_char != '\''))
    {
        /*  check the value only if it isn't N/A  */
        if ((strcmp(first_char, "n/a") != 0) && (strcmp(first_char, "N/A") != 0))
        {
            /*  we can't have multiple underscores in an unquoted value  */
            if (strstr(first_char, "__") != NULL)
            {
                (void)sprintf(message, "Multiple underscores in VALUE%s", element_prompt);
                status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
            }
    
            /*  an unquoted value cannot begin with an underscore  */
            if (*first_char == '_')
            {
                (void)sprintf(message, "First character is an underscore in VALUE%s", element_prompt);
                status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
            }
    
            /*  an unquoted value cannot end with an underscore  */
            if (*last_char == '_')
            {
                (void)sprintf(message, "Last character is an underscore in VALUE%s", element_prompt);
                status = OdlPrintMessage(message_fname,message_fptr, line_number, message, suppress_messages);
            }

            /*  the value may have a units expression  */
            if ((units_start = (char *) strchr(text, (int) '<')) != NULL)
	    {
                CopyString(save_units, units_start)
                *units_start = '\0';
		StripTrailing(text, ' ')
            }
            
            if (OdlDataType(text) == ODL_UNKNOWN)
            {
                (void)sprintf(message, "Unable to determine the data type of VALUE%s: \"%s\"", 
                                 element_prompt, first_char);
                status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
            }

            /*  validate the units expression, if any  */
            if (units_start != NULL)
            {
                /*  only one '<' char allowed in a units expression  */
                if (strchr(units_start+1, (int) '<') != NULL)
                {
                    (void)sprintf(message, "Embedded '<' character found in the UNITS expression: \"<%s\", for VALUE%s:  \"%s\"", 
                                     units_start+1, element_prompt, first_char);
                    status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
                }

                /*  find the closing char for the units expression  */
                units_end = (char *) strchr(units_start+1, (int) '>');

                /*  missing the closing '>' char in the units expression  */
                if (units_end == NULL)
                {
                    (void)sprintf(message, "Missing the closing '>' character in the UNITS expression: \"<%s\", for VALUE%s:  \"%s\"", 
                                     units_start+1, element_prompt, first_char);
                    status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
                }
                else
                /*  characters found after the closing '>' in the units exp  */
                if (units_end != last_char)
                {
                    (void)sprintf(message, "Extraneous characters found after the closing '>' character in the UNITS expression: \"<%s\", for VALUE%s:  \"%s\"", 
                                     units_start+1, element_prompt, first_char);
                    status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
                }
  
                /*  restore the value  */
                (void)strcat(text, " ");
                (void)strcat(text, save_units);
                LemmeGo(save_units)

            }  /*  End:  "if (units_start != NULL) ..."  */

        }  /*  End:  "if ((strcmp(first_char, "n/a") != 0) && ..."  */

    }  /*  End:  "if ((double_quote > ... else ... else ..."  */
    
    LemmeGo(message)

    return(status);

}  /*  End routine:  "OdlValidElement"  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlValidEndObjDesc (curr_object, equals, right_part, 
                             message_fname, message_fptr, line_number, group_type,
                             suppress_messages)
OBJDESC *curr_object;
char *equals;
char *right_part;
char *message_fname;
FILE *message_fptr;
long line_number;
unsigned short group_type;
int suppress_messages;

#else

short OdlValidEndObjDesc (OBJDESC *curr_object, char *equals, 
                                 char *right_part, char *message_fname, 
                                 FILE *message_fptr, long line_number, unsigned short group_type,
                                 int suppress_messages)

#endif
{
	short status = {TRUE};
	char errmsg[120];

	if (curr_object->parent == NULL)
	{
		/* Found an extra END_GROUP or END_OBJECT */
		(void)sprintf (errmsg,
		    "Encountered an extra END_%s - Ignored",
		    (group_type==ODL_OBJECT) ? "OBJECT" : "GROUP");
		status = OdlPrintMessage(message_fname, message_fptr, line_number, errmsg, suppress_messages);
	}

	if (*equals != '\0')
	{
		if (curr_object->is_a_group == group_type)
		{
			/* Make sure that the name in the END_OBJECT or END_GROUP part
                           matches the name that was given previously in the OBJECT
                           or GROUP part */

			if (*right_part != '\0')
			{
				if (strcmp(curr_object->class, right_part) != 0)
				{
					(void)sprintf (errmsg,
					    "END_%s = %s doesn't match %s = %s",
					    (group_type==ODL_OBJECT) ? "OBJECT" : "GROUP",
					    right_part,
					    (curr_object->is_a_group==ODL_OBJECT) ? "OBJECT" : "GROUP",
					    curr_object->class);
					status = OdlPrintMessage(message_fname, message_fptr, line_number, errmsg, suppress_messages);
				}
			}

			status = OdlValidIdentifier(right_part, "END_OBJECT class", 
			    message_fname, message_fptr,line_number, suppress_messages) && status;
		}
		else
		{
			/* Found an END_GROUP when expecting END_OBJECT, or vice versa */

			(void)sprintf (errmsg,
			    "Found END_%s when expecting END_%s - Ignored",
			    (group_type==ODL_OBJECT) ? "OBJECT" : "GROUP",
			    (curr_object->is_a_group==ODL_OBJECT) ? "OBJECT" : "GROUP");
			status = OdlPrintMessage(message_fname, message_fptr, line_number, errmsg, suppress_messages);
		}

	}  /*  End:  "if (*equals != '\0') ..."  */

	return(status);

}  /*  End routine:  "OdlValidEndObjDesc"  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlValidIdentifier (id_name, id_type, message_fname, message_fptr, line_number,
                          suppress_messages)

const char *id_name;
const char *id_type;
const char *message_fname;
FILE *message_fptr;
long line_number;
int suppress_messages;

#else

short OdlValidIdentifier (const char *id_name, const char *id_type, 
                          const char *message_fname, FILE *message_fptr, 
                          long line_number, int suppress_messages)

#endif
{
    char *message = NULL;
    char *c = {NULL};
    int i;
    short status = {TRUE};

    NewString(message, (TB_MAXLINE+(long)strlen(id_name)+(long)strlen(id_type)))

    if (id_name == NULL)
    {
        (void)sprintf(message, "%s identifier is missing", id_type);
        status = OdlPrintMessage(message_fname, message_fptr, line_number,message, suppress_messages);
    }
    else
    {
        char *local_id_name = strdup(id_name);
        if (!local_id_name)
            SayGoodbye();

        StripUnprintables(local_id_name)

        if (*local_id_name == '\0')
        {
            (void)sprintf(message, "%s identifier is missing", id_type);
            status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
        }
        else
        {
            if (! isalpha(*local_id_name))
            {
                (void)sprintf(message, 
                        "%s identifier:  \"%s\"  does not begin with a letter",
                        id_type, local_id_name);
                status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
            }
        
            for (c=local_id_name,i=0; *c != '\0'; ++c)
            {
                if ((*c != '_') && (! isalnum(*c))) ++i;
            }
        
            if (i > 0)
            {
                (void)sprintf(message, 
                        "%s identifier:  \"%s\"  contains %d embedded non-alphanumeric or \"_\" character", 
                        id_type, id_name, i);
                if (i > 1) (void)strcat(message, "s");
                status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
            }

        }  /*  End:  "if (*id_name == '\0') ... else ..."  */

        free(local_id_name);
    }  /*  End:  "if (id_name == NULL) ... else ..."  */
    
    LemmeGo(message)

    return(status);

}  /*  End routine:  "OdlValidIdentifier"  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlValidKwd (curr_object, left_part, equals, right_part, 
                   message_fname, message_fptr, line_number, suppress_messages)
OBJDESC *curr_object;
char *left_part;
char *equals;
char *right_part;
char *message_fname;
FILE *message_fptr;
long line_number;
int suppress_messages;

#else

short OdlValidKwd (OBJDESC *curr_object, char *left_part, char *equals,
                          char *right_part, char *message_fname, 
                          FILE *message_fptr, long line_number,
                          int suppress_messages)

#endif
{
    KEYWORD *keyword = {NULL};
    char *key = {NULL};
    char *message = NULL;
    short status = {TRUE};
    short sfdu_only = {FALSE};
    short found_keyword = {FALSE};

    NewString(message, (TB_MAXLINE+(long)strlen(left_part)+(long)strlen(right_part)))

    if (*left_part == '=')
    {
        *left_part = '\0';
        status = OdlPrintMessage(message_fname, message_fptr, line_number, 
                     "KEYWORD identifier is missing", suppress_messages);
    }
    else
    {
        if (*equals == '\0')
        {
            sfdu_only = ((strncmp(left_part, "NJPL", 4) == 0) ||
                            (strncmp(left_part, "CCSD", 4) == 0));

            if (! sfdu_only)
            {
                (void)sprintf(message, 
                    "Missing equals sign after KEYWORD identifier:  \"%s\"", left_part);
                status = OdlPrintMessage(message_fname, message_fptr, line_number, message, suppress_messages);
	    }
        }
        else
        {
            *equals = '\0';
            StripTrailing(left_part, ' ')
        }

        /*  ignore the first character if the keyword is a pointer  */
        key = (*left_part != '^') ? left_part : left_part + 1;

        status = OdlValidIdentifier(key, "KEYWORD", message_fname, 
                                    message_fptr, line_number, suppress_messages) && status;

        for (keyword=curr_object->first_keyword;
                ((keyword != NULL) && (! found_keyword));
                    keyword=keyword->right_sibling)
        {
            if (keyword->name != NULL)
                found_keyword = (strcmp(keyword->name, left_part) == 0);
        }

        if (found_keyword)
        {
            (void)sprintf(message, 
                    "Duplicate KEYWORD identifier:  \"%s\"", left_part);
            status = OdlPrintMessage(message_fname, message_fptr, line_number,message, suppress_messages);
        }
        
    }  /*  End:  "if (*left_part == '=') ... else ..."  */

    if (*right_part != '\0')
    {
        /*  what sort of value do we have?  */
        if ((*right_part != '{') && (*right_part != '('))
        {
            /*  we have a single element */
            status = OdlValidElement(right_part, message_fname, message_fptr, 
                                  line_number, (long)0, suppress_messages) && status;
        }
        else
        {
            /*  we have a value list  */
            status = OdlValidValueList(right_part, message_fname, message_fptr,
                                    line_number, suppress_messages) && status;
        }
    }
    else
    if (! sfdu_only)
    {
        (void)sprintf(message, 
                "KEYWORD identifier:  \"%s\"  is missing a VALUE",
                left_part);
        status = OdlPrintMessage(message_fname, message_fptr, line_number,message, suppress_messages);
    }

    LemmeGo(message)

    return(status);

}  /*  End routine:  "OdlValidKwd"  */


/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlValidObjDesc (curr_object, equals, right_part, 
                          message_fname, message_fptr, line_number, suppress_messages)

OBJDESC *curr_object;
char *equals;
char *right_part;
char *message_fname;
FILE *message_fptr;
long line_number;
int suppress_messages;

#else

short OdlValidObjDesc (OBJDESC *curr_object, char *equals, 
                              char *right_part, char *message_fname, 
                              FILE *message_fptr, long line_number,
                              int suppress_messages)

#endif
{
    short status = {TRUE};

    if (*equals == '\0')
    {
        status = OdlPrintMessage(message_fname, message_fptr, line_number, 
                     "Missing equals sign after OBJECT statement", suppress_messages);
    }

    StripTrailing(right_part, ' ')

    status = OdlValidIdentifier(right_part, "OBJECT class", 
                 message_fname, message_fptr,line_number, suppress_messages) && status;

    return(status);

}  /*  End routine:  "OdlValidObjDesc"  */



/*******************/
/*  Local Routine  */
/*******************/

#ifdef _NO_PROTO

short OdlValidValueList (text, message_fname, message_fptr, line_number, suppress_messages)

char *text;
char *message_fname;
FILE *message_fptr;
long line_number;
int suppress_messages;

#else

short OdlValidValueList (char *text, char *message_fname, 
                         FILE *message_fptr, long line_number,
                         int suppress_messages)

#endif
{
    char *first_char = {NULL};
    char *last_char = {NULL};
    char save_c;
    long i;
    short status = {TRUE};

    for (i=1,first_char=OdlValueStart(text); *first_char != '\0'; ++i)
    {
        /*  find the end of the current element  */
        last_char = OdlValueEnd(first_char);

        /*  save the next character and terminate the string  */
        save_c = *(++last_char);
        *last_char = '\0';

        /*  validate the current element  */
        (void)OdlValidElement(first_char,message_fname,message_fptr,line_number,i,
                              suppress_messages);

        /*  restore the character that was overwritten by the terminator  */
        *last_char = save_c;

        /*  find the start of the next element  */
        first_char = OdlValueStart(last_char);

    }  /*  End:  "for (i=1, ..."  */

    return(status);

}  /*  End routine:  "OdlValidValueList"  */




/*========================================================================*/
/*                                                                        */
/*                        Miscellaneous routines                          */
/*                                                                        */
/*========================================================================*/


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlWildCardCompare                                              */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine compares two text strings, one of which may have   */
/*      wildcard characters in it ('*').                                */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

unsigned short OdlWildCardCompare (wildcard_text, plain_text)

const char *wildcard_text;
const char *plain_text;

#else

unsigned short OdlWildCardCompare (const char *wildcard_text, 
                                   const char *plain_text)

#endif
{
    char *c = {NULL};
    char *substr = {NULL};
    char *tmp_str = {NULL};
    char *tmp_str2 = {NULL};
    char *text_start = {NULL};
    char save_it;
    unsigned long len; 
    unsigned short allrightythen = {FALSE};
    int  StartOfLine = 1;

    /*  see if we have anything to compare  */
    if ((wildcard_text != NULL) && (plain_text != NULL))
    {
        /*  all righty then, let's initialize some local variables  */
        allrightythen = TRUE;

        /*  copy the wildcard text  */
        CopyString(tmp_str, wildcard_text)

        /*  strip off leading and trailing quotes  */
        save_it = *tmp_str;
        if ((save_it == '\'') || (save_it == '"'))
        {
            StripLeading(tmp_str, save_it)
            StripTrailing(tmp_str, save_it)
            StripLeading(tmp_str, ' ')
            StripTrailing(tmp_str, ' ')
        }

        /*  copy the plain text  */
        CopyString(tmp_str2, plain_text)

        /*  strip off leading and trailing quotes  */
        save_it = *tmp_str2;
        if ((save_it == '\'') || (save_it == '"'))
        {
            StripLeading(tmp_str2, save_it)
            StripTrailing(tmp_str2, save_it)
            StripLeading(tmp_str2, ' ')
            StripTrailing(tmp_str2, ' ')
        }

        substr = tmp_str;
	/* if * is in the front of the substr, then don't check to see the the front of the substr and text_start match */
	if (*substr == '*')
        {
	    StartOfLine = 0;
        }
        text_start = tmp_str2;

        if (strchr(substr, '*') == NULL)
            allrightythen = (strcmp(text_start, substr) == 0);
        else
        {

            /*  we're going to break out the chunks of text between the */
            /*  wildcard caracters and try to find them one-by-one in   */
            /*  the plain text string.                                      */
            for(;;)
            {
                /*  locate the start of next substring  */
                for ( ; *substr == '*'; ++substr) ;
                if (*substr == '\0') break;

                /*  locate the end of the substring and save that address  */
                for (c=substr; ((*c != '*') && (*c != '\0')); ++c) ;
                save_it = *c;
                *c = '\0';
    
                /*  look for the substring in the un-wildcarded text  */
                if ((c = (char *)strstr(text_start, substr)) == NULL)
                {
                    allrightythen = FALSE;
                    break;
                }
		else if (StartOfLine)
		{
		    if (c != text_start)
		    {
			allrightythen = FALSE;
			break;
		    }
		    StartOfLine = 0;
		}

                /*  prepare for the next search  */
                len = strlen(substr);
                substr += len;
                *substr = save_it;
                text_start = c + len;
                if ((*substr == '\0') && (*text_start != '\0'))
                {
                    allrightythen = FALSE;
                    break;
                }

            }

	}  /*  End:  "if (strchr(substr, '*') == NULL) ... else ..."  */

        LemmeGo(tmp_str)
        LemmeGo(tmp_str2)

    }  /*  End:  "if ((wildcard_text != NULL) && ..."  */

    return(allrightythen);

}  /*  End:  "OdlWildCardCompare"  */




/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlTraverseTree                                                 */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the next object in an ODL tree, stopping   */
/*      when it has traversed the entire tree as defined by the         */
/*      root_level parameter.                                           */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

OBJDESC *OdlTraverseTree (curr_object, root_level)

OBJDESC *curr_object;
unsigned long root_level;

#else

OBJDESC *OdlTraverseTree (OBJDESC *curr_object, unsigned long root_level)

#endif
{
    OBJDESC *obj = {NULL};
    OBJDESC *next_object = {NULL};

    if (curr_object != NULL)
    {
        /*  start search with current object's children  */
        if (curr_object->first_child != NULL)
            next_object = curr_object->first_child;
        else
        /*  start search with current object's right sibling  */
        if (curr_object->right_sibling != NULL)
            next_object = curr_object->right_sibling;
        else
        /*  move up parent list until we find one with a right sibling  */
        {
            for (next_object=NULL,obj=curr_object->parent; 
                    (obj != NULL); obj=obj->parent)
            {
                if (obj->level <= root_level) 
                    break;
                else
                if (obj->right_sibling != NULL)
                {
                    next_object = obj->right_sibling;
                    break;    
                }
            }

        }  /*  End:  "if (curr_object->first_child ... else ..."  */

    }  /*  End:  "if (curr_object != NULL) ..."  */

    return(next_object);

}  /*  End routine:  "OdlTraverseTree"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlFirstWord                                                    */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a pointer to the first word in a string    */
/*      of text.  A word is anything that begins with a printable       */
/*      ASCII character.                                                */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlFirstWord(text)

char *text;

#else

char *OdlFirstWord(char *text)

#endif
{
    char *c = {NULL};
    for (c=text; 
          ((c != NULL) && ((*c <= ' ') || (*c > '~')) && (*c != '\0')); ++c) ;
    return(c);

}  /*  End routine:  "OdlFirstWord"  */


	
/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlNextWord                                                     */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the next word in a string of text.  It     */
/*      locates the end of the current word and skips over whitespace   */
/*      to find the start of the next.  A word is anything that begins  */
/*      with a printable ASCII character.                               */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlNextWord(text)

char *text;

#else

char *OdlNextWord( char *text)

#endif
{
    char *c = {NULL};
    for (c=text; 
          ((c != NULL) && (*c > ' ') && (*c <= '~') && (*c != '\0')); ++c) ;
    return(OdlFirstWord(c));

}  /*  End routine:  "OdlNextWord"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlValueStart                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the start of a value within a set or       */
/*      sequence.  A value is anything that doesn't begin with a brace, */
/*      comma, paren, or blank.                                         */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlValueStart(text)

char *text;

#else

char *OdlValueStart(char *text)

#endif
{
    /*  find a character that is not a brace, paren, comma, or blank  */
    return(text + strspn(text, "{}(), \n"));

}  /*  End routine:  "OdlValueStart"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlValueEnd                                                     */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine finds the end of a value within a set or sequence. */
/*      A value is anything that doesn't begin with a brace, comma,     */
/*      paren, or blank.  It finds the end by locating a comma, paren,  */
/*      or brace, and backing up over trailing whitespace.              */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlValueEnd(text)

char *text;

#else

char *OdlValueEnd( char *text)

#endif
{
    char *c = {NULL};

    /*  find a character that is a brace, paren, or comma  */
    c = strpbrk(text, "{}(),");

    /*  backup over any trailing blanks  */
    for (--c; ((c > text) && ((*c == ' ') || (*c == '\0'))); --c) ;

    return(c);

}  /*  End routine:  "OdlValueEnd"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlValueRowStart                                                */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the start of a set or sequence within a    */
/*      set or sequence.  A set begins with a brace and a sequence      */
/*      begins with a paren.                                            */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlValueRowStart(text)
char *text;

#else

char *OdlValueRowStart( char *text)

#endif
{
    /*  find a character that is not a brace or paren  */
    return(text + strspn(text, "{}()\n"));

}  /*  End routine:  "OdlValueRowStart"  */


/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlValueRowEnd                                                  */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine locates the end of a set or sequence within a      */
/*      set or sequence.  A set begins with a brace and a sequence      */
/*      begins with a paren.  It finds a closing brace or paren and     */
/*      backs up over any whitespace.                                   */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlValueRowEnd(text)
char *text;	

#else

char *OdlValueRowEnd( char *text)

#endif
{
    char *c = {NULL};

    /*  find a character that is a brace or paren  */
    c = strpbrk(text, "{}()\n");

    /*  backup over any trailing blanks or commas  */
    for (--c; ((c > text) && ((*c == ' ') || (*c == ',') || (*c == '\0'))); --c) ;

    return(c);

}  /*  End routine:  "OdlValueRowEnd"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlDataType                                                     */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*      04-09-97    Added Landsat 7 Date and Time 	                */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine determines the ODL data type of a text string.     */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

unsigned short OdlDataType (text)

char *text;

#else

unsigned short OdlDataType (char *text)

#endif
{
    char *c = {NULL};
    char *u = {NULL};
    char *tempstr = {NULL};
    char *last_pound = {NULL};
    int base;
    unsigned long decimal_count = {0};
    unsigned long hyphen_count = {0};
    unsigned long colon_count = {0};
    unsigned long t_count = {0};
    unsigned short has_sign = {FALSE};
    unsigned short type = {ODL_UNKNOWN};
    unsigned short exponent_type = {ODL_UNKNOWN};
    static const char *valid_chars[17] = {
        "",   /* base 0 */
        "0",   /* base 1 */
        "01",   /* base 2 */
        "012",   /* base 3 */
        "0123",   /* base 4 */
        "01234",   /* base 5 */
        "012345",   /* base 6 */
        "0123456",   /* base 7 */
        "01234567",   /* base 8 */
        "012345678",   /* base 9 */
        "0123456789",   /* base 10 */
        "0123456789A",   /* base 11 */
        "0123456789AB",   /* base 12 */
        "0123456789ABC",   /* base 13 */
        "0123456789ABCD",   /* base 14 */
        "0123456789ABCDE",   /* base 15 */
        "0123456789ABCDEF"};  /* base 16 */

    if (text != NULL)
    {
        /*  make a copy of the string  */
        CopyString(tempstr, text)

        /*  remove any units  */
        if ((u = (char *)strchr(tempstr, '<')) != NULL) *u = '\0';

        /*  clean up blanks and convert to upper case  */
        StripLeading(tempstr, ' ')
        StripTrailing(tempstr, ' ')
        UpperCase(tempstr)

        c = tempstr;

        type = ODL_SYMBOL;

        /*  See if we have a double quoted text value  */
        if (*c == '"')
            type = ODL_TEXT;
        else
        /*  See if we have a single quoted text value  */
        if (*c == '\'')
            type = ODL_SYMBOL;
        else
        /*  See if we have a sequence of values  */
        if (*c == '{')
            type = ODL_SEQUENCE;
        else
        /*  See if we have a set of values  */
        if (*c == '(')
            type = ODL_SET;
        else
        /*  See if we have any embedded blanks  */
        if (strchr(c, ' ') != NULL)
            type = ODL_UNKNOWN;
        else
        /*  See if we have an integer, real, date, or date-time value  */
        {
            /*  we're going to check each character  */
            for ( ; *c != '\0'; ++c)
            {
                /*  may have a number  */
                if (isdigit(*c))
                {
                    if (c == tempstr)
                        type = ODL_INTEGER;
                    else
                    if (has_sign && (c == (tempstr + 1)))
                        type = ODL_INTEGER;
                }
                else
                /*  we may have a real number or a date-time  */
                if (*c == '.')
                {
                    /*  we may have a real number  */
                    if (type == ODL_INTEGER)
                        type = ODL_REAL;
                    else
                    /*  date-times can only have one decimal point  */
                    if ( (type == ODL_DATE_TIME) || (type == ODL_L7_DATE_TIME_FRAC) )
                    {
                        if ((++decimal_count) > 1)
                            type = ODL_UNKNOWN;
		    }
                    else
		    /* Landsat 7 date and time fraction */
		    if (type == ODL_L7_DATE_TIME)
		       type = ODL_L7_DATE_TIME_FRAC;
                    else
                    /*  we may have a real number  */
                    if (c == tempstr)
                        type = ODL_REAL;
                    else
                    /*  we may have a signed real number  */
                    if (has_sign && (c == (tempstr + 1)))
                        type = ODL_REAL;
                    else
                        type = ODL_UNKNOWN;
                }
                else
                /*  we may have a real number in scientific notation */
                if (*c == 'E')
                {
                    /*  only valid if we thought we had an int or real  */
                    if ((type == ODL_INTEGER) || (type == ODL_REAL))
                    {
                        /*  check out the exponent  */
                        exponent_type = (unsigned short) OdlDataType((c+1));

                        /*  we have a real number only if the exponent  */
                        /*  is real or int                              */
                        if ((exponent_type == ODL_REAL) || (exponent_type == ODL_INTEGER))
                            type = ODL_REAL;
                        else
                            type = ODL_UNKNOWN;

                        break;
                    }
                    else
                        if (type != ODL_SYMBOL)
                            type = ODL_UNKNOWN;
                }
                else
                /*  we may have a signed number  */
                if (*c == '+')
                {
                    /*  this had better be the first character  */
                    if (c != tempstr)
                        type = ODL_UNKNOWN;
                    else
                        has_sign = TRUE;
                }
                else
                /*  we may have a date or a signed number */
                if (*c == '-')
                {
                    /*  this had better be the first character  */
                    if (c == tempstr)
                        has_sign = TRUE;
                    else
                    /*  a date can have at most two hyphens  */
                    if ((++hyphen_count) > 2)
                        type = ODL_UNKNOWN;
                    else
                    /*  we thought we had an integer  */
                    if (type == ODL_INTEGER)
                        type = ODL_DATE;
                    else
                    /*  if it wasn't an int and it wasn't a date ...  */
                    if (type != ODL_DATE)
                        type = ODL_UNKNOWN;
                }
                else
                /*  we may have a date-time  */
                if (*c == 'T')
                {
                    /*  we thought we had a date  */
                    if (type == ODL_DATE)
                        type = ODL_DATE_TIME;
                    else
                    /*  a date-time may only have one 'T'  */
                    if (type == ODL_DATE_TIME)
                    {
                        if ((++t_count) > 1)
                            type = ODL_UNKNOWN;
		    }
                    else
                    /*  must be a symbol  */
                    if (type != ODL_SYMBOL)
                        type = ODL_UNKNOWN;
                }
                else
                /*  we may have a date-time  */
                if (*c == 'Z')
                {
                    /*  only a date time may contain a 'Z'  */
                    if (type == ODL_DATE_TIME)
                    {
                        /*  it had better be the last char in the string  */
                        if (*(c + 1) != '\0')
                            type = ODL_UNKNOWN;
		    }
                    else
                    if (type != ODL_SYMBOL)
                        type = ODL_UNKNOWN;
                }
                else
                /*  we may have a date-time  */
                if (*c == ':')
                {
		    /* we have a Landsat 7 date and time */
		    if (type == ODL_INTEGER)
		        type = ODL_L7_DATE_TIME;
		    /*  only a date-time may contain colons  */
		    if ( (type != ODL_DATE_TIME) && (type != ODL_L7_DATE_TIME) )
		        type = ODL_UNKNOWN;
                    else
                    /*  there can't be more than two of them in date  */
                    if ( ((++colon_count) > 2) && (type != ODL_L7_DATE_TIME ) )
                        type = ODL_UNKNOWN;
                    else
		    if ( ( colon_count > 4 ) && (type == ODL_L7_DATE_TIME) )
		       type = ODL_UNKNOWN;
                    else
                    /*  characters on either side must be digits  */
                    if ((! isdigit(*(c-1))) || (! isdigit(*(c+1))))
                        type = ODL_UNKNOWN;
                    else
                    /*  decimal points can't occur before the last colon  */
                    if (decimal_count != 0)
                        type = ODL_UNKNOWN;
                }
                else
                /*  we may have a non-decimal integer  */
                if (*c == '#')
                {
                    /*  we didn't think it WAS an integer  */
                    if (type != ODL_INTEGER)
                        type = ODL_UNKNOWN;
                    else
                    /*  the base can't be signed  */
                    if (has_sign)
                        type = ODL_UNKNOWN;
                    else
                    /*  missing the closing '#' character  */
                    if ((last_pound = (char *)strchr(c+1, '#')) == NULL) 
			 type = ODL_UNKNOWN;
		    else
                    /*  closing '#' char is not at the end  */
                    if (*(last_pound + 1) != '\0')
                        type = ODL_UNKNOWN;
                    else
                    {
                        /*  looks good so far, but we have to make sure  */
                        /*  stuff between the '#' characters is valid    */
                        /*  with the specified base                      */

                        /*  isolate the base  */
                        *c = '\0'; 
                        base = atoi(tempstr);

                        /*  ignore the integer's sign  */
                        ++c;
                        if ((*c == '+') || (*c == '-')) ++c;

                        /*  isolate the number part  */
                        *last_pound = '\0';

                        /*  valid bases are 2 through 16, inclusive  */
                        if ((base < 2) || (base > 16))
                            type = ODL_UNKNOWN;
                        else
                        /*  look for invalid digits for the specified base  */
                        if (c[strspn(c, valid_chars[base])] != '\0')
                            type = ODL_UNKNOWN;
                    }

                    if (type != ODL_UNKNOWN) 
                    {
                        type = ODL_INTEGER;
                        break;
		    }
                }
                else
                /*  we may have an unquoted symbol  */
                if (isalpha(*c) || (*c == '_'))
                {
                    if  (type != ODL_SYMBOL)
                        type = ODL_UNKNOWN;
                }
                else
                /*  we havn't got a clue  */
                    type = ODL_UNKNOWN;

                if (type == ODL_UNKNOWN) break;

            }  /*  End:  "for ( ; ..."  */

            if (has_sign && (type != ODL_INTEGER) && (type != ODL_REAL))
                type = ODL_UNKNOWN;

        }  /*  End:  "if (*c == '(') ... else ..."  */

        LemmeGo(tempstr)

    }  /*  End:  "if (text != NULL) ..."  */

    return(type);

}  /*  End:  "OdlDataType"  */

/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlTypeString                                                   */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine determines the ODL data type of a text string.     */
/*                                                                      */
/************************************************************************/

#ifdef _NO_PROTO

char *OdlTypeString (type, type_string)

unsigned short type;
char *type_string;

#else

char *OdlTypeString (unsigned short type, char *type_string)

#endif
{
    static char local_type_string [TB_MAXLINE];

    switch (type)
    {
        case ODL_INTEGER   : (void)strcpy(local_type_string, "INTEGER");
                             break;

        case ODL_REAL      : (void)strcpy(local_type_string, "REAL");
                             break;

        case ODL_SYMBOL    : (void)strcpy(local_type_string, "SYMBOL");
                             break;

        case ODL_TEXT      : (void)strcpy(local_type_string, "TEXT");
                             break;

        case ODL_DATE      : (void)strcpy(local_type_string, "DATE");
                             break;

        case ODL_DATE_TIME : (void)strcpy(local_type_string, "DATE-TIME");
                             break;

        case ODL_SEQUENCE  : (void)strcpy(local_type_string, "SEQUENCE");
                             break;

        case ODL_SET       : (void)strcpy(local_type_string, "SET");
                             break;

        default            : (void)strcpy(local_type_string, "UNKNOWN");
                             break;

    }  /*  End:  "switch (type) ..."  */

    if (type_string != NULL) (void)strcpy(type_string, local_type_string);

    return(local_type_string);

}  /*  End:  "OdlTypeString"  */



/************************************************************************/
/*                                                                      */
/*  Component:                                                          */
/*                                                                      */
/*      OdlTempFname                                                    */
/*                                                                      */
/*  Author:                                                             */
/*                                                                      */
/*      David P. Bernath (Jet Propulsion Laboratory)                    */
/*                                                                      */
/*  Version:                                                            */
/*                                                                      */
/*      1.0    March 31, 1994                                           */
/*                                                                      */
/*  Change History:                                                     */
/*                                                                      */
/*      03-31-94    Original code                                       */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*                                                                      */
/*  Description:                                                        */
/*                                                                      */
/*      This routine returns a unique file name, depending on the       */
/*      system on which the code was compiled.                          */
/*                                                                      */
/************************************************************************/

char *OdlTempFname()
{
    FILE *fptr = {NULL};
    char *fname = {NULL};
    char temp_str  [TB_MAXPATH + TB_MAXFNAME];
    char base_name [TB_MAXPATH + TB_MAXFNAME];

    (void)strcpy(base_name, "tmp.tmp");

#ifdef SUN_UNIX
    (void)tmpnam(temp_str);
    (void)strcpy( base_name, temp_str);  /* Bug fix 11/2/94 SM                     */
                                   /* Was:    (void)sprintf(base_name, "~/%s.tmp", */
                                   /*                 temp_str);             */
#endif

#if (defined( VAX) || defined( ALPHA_VMS))
    (void)tmpnam(temp_str);
    (void)sprintf(base_name, "sys$login:%s.tmp", temp_str);
#endif

#ifdef MAC_THINK
    (void)tmpnam(temp_str);
    (void)sprintf(base_name, "%s.tmp", temp_str);
#endif

#ifdef MSDOS
    {
        time_t t;
        t = (time_t) time(NULL);
        (void)sprintf(base_name, "C:\\%ld", t);
        base_name[8] = EOS;
        (void)strcat(base_name, ".tmp");
    }
#endif

    CopyString(fname, base_name)

    if ((fptr = (FILE *) fopen(fname, "w")) == NULL)
        LemmeGo(fname)
    else
        CloseMe(fptr)

    return(fname);

}  /*  End:  "OdlTempFname"  */


short CheckBalance(text)
char *text;
{
    long quote_nesting = 0;
    long brace_nesting = 0;
    long paren_nesting = 0;
    char *c = {NULL};
    char *c1 = {NULL};

    c1 = (char *) strchr(text, '=');
    c = OdlFirstWord(c1 + 1);

    if ((*c == '(') || (*c == '{'))
        (void)OdlNestingLevel(c,&brace_nesting,&paren_nesting);
    else
        if (*c == '"')
        {
            for (; *c != '\0'; ++c)
            {
                if (*c == '"') 
                {
                    if (quote_nesting == 0)
                        quote_nesting = 1;
                    else
                        quote_nesting = 0;
	        }
            }
        }

    return((brace_nesting + paren_nesting + quote_nesting) == 0);
}

/*========================================================================*/
/*                                                                        */
/*                       End of lablib 3.0 stuff                          */
/*                                                                        */
/*========================================================================*/

