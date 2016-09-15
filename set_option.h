#ifndef SET_OPTION_H
#define SET_OPTION_H






void set_options(int argc, char **argv, int *opts);
void init_config(int *opts, CONF_t *init);
void Root_Bcat_data(CONF_t *init, int root);

#endif
