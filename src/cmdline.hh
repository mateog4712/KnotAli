#ifndef CMDLINE_H
#define CMDLINE_H
#include <string>


/** @brief Where the command line options are stored */
struct args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  const char *verbose_help; /**< @brief Turn on verbose output help description.  */
  const char *input_type_help; /**< @brief Give input file type help description.  */
  const char *output_file_help; /**< @brief Give output file help description.  */
  const char *stacking_help; /**< @brief Give stacking help description.  */
  const char *threads_help; /**< @brief Give threads help description.  */
  const char *pseudoknot_help; /**< @brief Give pseudoknot help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int input_type_given ;	/**< @brief Whether input file was given.  */
  unsigned int output_file_given ;	/**< @brief Whether output file was given.  */
  unsigned int stacking_given ;	/**< @brief Whether stacking was given.  */
  unsigned int threads_given ;	/**< @brief Whether threads was given.  */
  unsigned int pseudoknot_given ; /**< @brief Whether pseudoknot was given.  */

  char **inputs ; /**< @brief unnamed options (options without names) */
  unsigned inputs_num ; /**< @brief unnamed options number */
} ;

extern std::string type;
extern std::string file;
extern int numThreads;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];


/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,struct args_info *args_info);


/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv, struct args_info *args_info);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct args_info *args_info);

/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct args_info *args_info);






#endif