/*ms----------------------------------------------------------------------------

   tilesig.c

   Purpose:
      To calculate approximate original sigma nought value for all subtiles 
      in a RAMS final tile and write them to a full 32 bit output mosaic.

   Procedures:
      Depend   - To check program dependencies and collect input data
      usage     - To print a usage message 
      GetSigma0 - To get value from line & sample that a selected x,y map to.
      ParseArgs - To set defaults and parse the command line.

   Description:
      This program takes each pixel value from a final image subtiles and
      reverses all of the corrections applied by the RAMS software after the
      orthorectification, thereby achieving an approximate sigma nought value.

      The source was ripped out of Petes GETSIG8.C found in the cdrom directory.

      The program must reverse the effects of the following:
      grand_rad - Block to block radiometric adjustment
      grand_geo - Block to block geometric adjustment
      rad_bal   - Radiometric adjustment between frames in block

   Interface: tilesig -out output_image
         [-h]      - (help) print usage
         [-db]      - print additional internal info

   Environment
      This program needs to have access to the following files:

        MASTER.TXT
        BLOCKS.KEY
        FRAMES.KEY

   Author:   ripped apart and recombined by millerjd

   Algorithm:
      This program gets information about the source of each pixel from the
      index image INDICES.DIR/<name>.IDX.  This provides block and frame
      information.  The grand radiometric adjustment is then reversed using 
      the radiometric equations from the BLOCKS.KEY file.  Before applying 
      the frame radiometric adjustment from the FRAMES.KEY file, the (x,y) 
      location of the pixel prior to the grand geometric adjustment must be 
      found by inverting the geometric equation from the BLOCKS.KEY file.

----------------------------------------------------------------------------me*/

/* ---- Include Files ---- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>

#define MAMM 1

#define INDEX_DIR "IMGINDEX.DIR"

#define NO_DATA_VAL  -9999
#define OUT_NULL     -32767
#define DATA_SCALE    1638.35
#define OFFSET        30.0
#define OUT_OFFSET    32766

/* ---- Local data types ---- */

typedef struct {           /* command line options...      */
   char   *output_file;     /* output file name            */
   char   *index_file;      /* index file name            */
   double  map_x;           /* user specified coordinate  */
   double  map_y;
   int     do_point;        /* flag that user gave point   */
   int     depend;          /* Was "-depend" specified?    */
   int     debug;           /* Was "-db" specified?        */
   int     help;            /* Was "-h" specified?         */
} Options_t;

typedef struct IntXY_t {
   int x, y;
} IntXY_t;

typedef struct DoubleXY_t {
   double x, y;
} DoubleXY_t;

typedef struct coeff_s {
   double value;
   struct coeff_s *next;
} coeff_t;

typedef struct {
   coeff_t *firstcoeff;
   int     n_coeffs;
} coeffs_t;

typedef struct {
    DoubleXY_t  map_xy;     /* map coords of chip center in parent image */
    DoubleXY_t  ul;         /* coords of upper-left chip corner in parent */
    double      avg;        /* chip magnitude */
    double      target;     /* target avg for balancing */
} EdgeTie_t;

typedef struct {
    EdgeTie_t   *tie;       /* array of edge ties */
    int         n_ties;     /* number of ties in array */
    IntXY_t     size;       /* size of chip */
    double      spacing;    /* spacing between ties */
} EdgeTies_t;


typedef struct {           /* Block definition */
   int           id;              /* block id associated with frame */
   char         *name;            /* block name from master file    */
   coeffs_t      blk_offset;      /* list of coefficients           */
   coeffs_t      blk_scale;       /* list of coefficients           */
   coeffs_t      blk_geom;        /* list of coefficients           */
   EdgeTies_t    blk_edgeties;    /* edge balancing pts for block   */
   } Block_t;

typedef struct {           /* Frame definition            */
   int           index;            /* index number from index file   */
   char         *name;             /* frame name from master file    */
   int           block_id;         /* id of associated block         */
   double        min_pwr;          /* conversion min power parameter */
   double        cnvt_scale;       /* conversion scale paramter      */
   coeffs_t      frm_offset;       /* list of coefficients           */
   coeffs_t      frm_scale;        /* list of coefficients           */
   EdgeTies_t    frm_edgeties;     /* edge balancing pts for frame   */
   Block_t      *block;            /* block reference for index      */
   } Frame_t;

typedef struct {           /* Subtile definition            */
   char  name[12];          /* base subtile name             */
   double min_x;           /* map extents                   */
   double max_x;
   double min_y;
   double max_y;
   int    img_ul_x;        /* pixel extents in output file  */
   int    img_ul_y;
   int    img_lr_x;
   int    img_lr_y;
   int    index_ul_x;        /* pixel extents in output file  */
   int    index_ul_y;
   short         *buf;
   unsigned char *i_buf; 
   }  Subtile_t;

typedef struct {           /* Subtile definition            */
   char  *name;            /* base subtile name             */
   float *buf;
   double min_x;           /* map extents                   */
   double max_x;
   double min_y;
   double max_y;
   int    size_x;
   int    size_y;
   int    fd;
   }  Image_t;

typedef struct {           /* data to be passed around like some cheap tart */
 /*  Variables to be changed on the fly */
   double      x;               /* current x coord to convert         */
   double      y;               /* current y coord to convert         */
   double      geo_x;           /* geometric map x adjustment         */
   double      geo_y;           /* geometric map y adjustment         */
   double      s0;              /* returned sigma nought value        */
   double      value;           /* input image value                  */
   int         index_value;     /* image index value                  */
 /* this stuff will stay put */
   double      image_res;       /* image pixel spacing                */
   double      index_res;       /* index image pixel spacing          */
   double      tile_size;       /* edge length of square subtile      */
   int         image_size;      /* number of pixels/lines on a side   */
   int         index_size;      /* number of pixels/lines on a side   */
   Frame_t    *frames;          /* array of frames                    */
   int         n_frames;        /* number of frames expected in index */
   Block_t    *blocks;          /* array of blocks                    */
   int         n_blocks;        /* number of blocks in tile           */
   Subtile_t  *subs;            /* array of subs to be calculated     */
   int         n_subs;          /* number of subtiles in tile         */
   Image_t    *output_image;    /* output data */
   Image_t    *index_image;     /* output data */
} Data_t;

/* ---- Function Prototypes ---- */

/* user interface */
void usage(char *cmd);
int ParseArgs(int argc, char *argv[], Options_t *options);

/* upfront collection of parameters */

int Depend( Options_t *options, Data_t *data);
coeffs_t *get_coeffs(char *line, Options_t *);
int calculate_output_parameters(Data_t *data, Options_t *options);

/* bulk of the work goes here */
int GetSigma0 (Options_t *options, Data_t *data);
int calculate_sub(Subtile_t *sub, Options_t *options, Data_t *data);
double ApplyEqnAtPt(double xx, double yy, coeffs_t *coeffs);
double FindOffsetAtPt(double xx, double yy, EdgeTie_t *tie_array, int n_ties, double spacing);

/* write to output */
int prepare_output(Data_t *data, Options_t *options);
int write_sub(Subtile_t *sub, Data_t *data, Options_t *options);
int give_head(Data_t *data, Options_t *options);

/*fs----------------------------------------------------------------------------

    Procedure:   main

    Purpose:   Get sigma nought value for an entire tile

    Exits:   Exit status is 0 on success, 1 on failure

----------------------------------------------------------------------------fe*/

int main(int argc, char *argv[])
{
   Options_t *options;
   Data_t    *data;
   int        ii;

   options = (Options_t *)calloc(1, sizeof(Options_t));
   data    = (Data_t *)calloc(1, sizeof(Data_t));

   /* ---- Parse command line ---- */

   if (ParseArgs(argc, argv, options)) {
      printf("%s: error parsing argument string\n", argv[0]);
      usage(argv[0]);
      exit (1);
      }


   /* ---- collect dependencies: ---- */

   if (Depend(options, data)) {
      printf("%s: error collecting parameters\n", argv[0]);
      exit (1);
      }

   if(options->depend)
      exit(0);

   if(options->do_point == 0)
      prepare_output(data, options);

   for(ii = 0; ii < data->n_subs; ii++) {
      calculate_sub(&data->subs[ii], options, data);
/*
      if(calculate_sub(&data->subs[ii], options, data) ){
         printf("%s: error generating output for %s\n", argv[0], data->subs[ii].name);
         exit (1);
         }
*/
      }

   if(options->do_point == 1)
       exit(0);

   close(data->output_image->fd);
   if(data->index_image)
      close(data->index_image->fd);

   if(give_head(data, options) ){
      printf("%s: error writing image header", argv[0]);
      exit (1);
      }

   /* ---- Return success ---- */

   exit(0);
}


/*fs----------------------------------------------------------------------------

    Procedure:   ParseArgs

    Purpose:   To set parse the command line.

    Arguments:   int     argc      - # of arguments on command line.
      char     *argv[]  - array of command line arguments.
      Options_t *options - Pointer to command line options.

    Returns:   Returns options by reference.

      Returns 0 on success. Otherwise, returns 1 on failure.

----------------------------------------------------------------------------fe*/

int ParseArgs(int argc, char *argv[], Options_t *options)
{
   int ii;

   if(argc < 2) {
      printf("insufficient arguments\n");
      return 1;
      }

   for(ii = 1; ii < argc; ii++) {
      if(!strcmp(argv[ii], "-out")) {
         ii++;
         options->output_file = argv[ii];
         continue;
         }
      if(!strcmp(argv[ii], "-index")) {
         ii++;
         options->index_file = argv[ii];
         continue;
         }
      if(!strcmp(argv[ii], "-point")) {
         ii++;
         sscanf( argv[ii], "%lf", &options->map_x);
         ii++;
         sscanf( argv[ii], "%lf", &options->map_y);
         options->do_point = 1;
         continue;
         }
      if(!strcmp(argv[ii], "-db")) {
         ii++;
         sscanf(argv[ii], "%d", &options->debug);
         continue;
         }
      if(!strcmp(argv[ii], "-depend")) {
         options->depend = 1;
         }
      if(!strcmp(argv[ii], "-h")) {
         usage(argv[0]);
         exit(0);
         }
     }

   if(options->do_point)
      return 0;

   if(options->output_file == NULL) {
      printf("output file not specified\n");
      usage(argv[0]);
      exit(1);
      }

   return 0;
}

/*fs----------------------------------------------------------------------------

    Procedure:   usage

    Purpose:   print usage 

    Arguments:   char *cmd    - name of command 

    Returns:   No return value.

----------------------------------------------------------------------------fe*/

void usage(char *cmd)
{
   printf("\n%s: generates original sigma nought values for an entire tile.\n\n", cmd);
   printf("     Optionally it will generate a corresponding index mosaic.\n");
   printf("     If given a single point, it will print the sigma_0 for\n");
   printf("     that point and no output file will be generated\n\n");
   printf( "  %s [-out output_file]\n\n", cmd);
   printf( "    -index index_file\n");
   printf( "    -point <x y>         - calculate for sigma_0 at given point\n");
   printf( "    -db    <print_level> - set debug output level\n");
   printf( "    -h                   - print usage\n\n");
}

/*fs----------------------------------------------------------------------------

    Procedure:   int Depend(options, data)

    Purpose:   To read program dependencies

      A switch "-depend" could be used to check the dependencies
      of the program and print the paths of key files known to
      be missing.  Currently, any missing dependencies result in
      an error.

    Returns:   Returns 0 on success or 1 on failure.

----------------------------------------------------------------------------fe*/
int Depend(Options_t *options, Data_t *data)
{
   //static char fn[] = "Depend";
   char name[13];
   char path[256], line[512];
   char *strptr;
   int ii, jj;
   int  do_ties = 0;
   FILE *fp;

   Subtile_t *sub;
   Frame_t   *frame;
   Block_t   *block;
   EdgeTie_t *tie;

   if(options->debug > 0)
       printf("collecting dependencies\n");

   /* ---- Init ---- */

   fp = NULL;

  /* Read MASTER.TXT file */

   if(options->debug > 5)
       printf("reading MASTER.TXT\n");

   fp = fopen("MASTER.TXT", "r");
   if(fp == NULL) {
      printf("missing MASTER.TXT\n");
      return 1;
      }

   ii = 0;
   while(fgets(line, 511, fp)) {
      if(strstr(line, "Image pixel spacing")) {
          strptr = strstr(line, ":");
          strptr++;
          sscanf(strptr, "%lf", &data->image_res);
          if(options->debug >= 10)
              printf("Image pixel spacing = %lf\n", data->image_res);
          }
      if(strstr(line, "Index tile pixel spacing")) {
          strptr = strstr(line, ":");
          strptr++;
          sscanf(strptr, "%lf", &data->index_res);
          if(options->debug >= 10)
              printf("Index tile pixel spacing = %lf\n", data->index_res);
          }
      if(strstr(line, "Subtile size")) {
          strptr = strstr(line, ":");
          strptr++;
          sscanf(strptr, "%lf", &data->tile_size);
          if(options->debug >= 10)
              printf("Subtile size = %lf\n", data->tile_size);
          data->image_size = data->tile_size / data->image_res;
          data->index_size = data->tile_size / data->index_res;

          if(options->debug >= 10) {
             printf("image size = %d\n", data->image_size);
             printf("index size = %d\n", data->index_size);
             }
          }

      if(strstr(line, "Number of sub-tiles")) {
          strptr = strstr(line, ":");
          strptr++;
          sscanf(strptr, "%d", &data->n_subs);
          if(data->n_subs == 0) {
             printf("No subitles in MASTER.TXT\n");
             return 1;
             }
          data->subs = (Subtile_t *)calloc(data->n_subs, sizeof(Subtile_t));
          if(options->debug >= 10) 
             printf("number of subs = %d\n", data->n_subs);
          }

      if(strstr(line, "Subtile") && data->subs) {
           strptr = strstr(line, ":");
           strptr++;
           sub = &data->subs[ii];
           sscanf(strptr, "%s %lf %lf %lf %lf", sub->name, &sub->min_x,
                  &sub->min_y, &sub->max_x, &sub->max_y);
           if(options->debug >= 15) {
              printf("%d: %s %.0lf %.0lf %.0lf %.0lf\n", ii, sub->name, sub->min_x,
                  sub->min_y, sub->max_x, sub->max_y);
              }
           ii++;
           }
       }

   fclose(fp);

  /* Read Frames.key file */

   sprintf(path, "%s/FRAMES.KEY", INDEX_DIR);
   fp = fopen(path, "r");
   if(fp == NULL) {
      printf("missing %s\n", path);
      return 1;
      }

   if(options->debug >= 5)
       printf("reading %s\n", path);

 /* count up the number of index frame references and allocate array of frames */

   ii = 0;
   while(fgets(line, 511, fp)) {
      if(strstr(line, "Frame Index")) ii++;
      }
   rewind(fp);
   data->frames = (Frame_t *)calloc(ii, sizeof(Frame_t));
   data->n_frames = ii;
   ii = 0;

   if(options->debug >= 5)
      printf("initializing %d index entries\n", data->n_frames);

   while(fgets(line, 511, fp)) {

      if(strstr(line, "Frame Index")) { // ready to go with new frame
         frame = &data->frames[ii];
         strptr = strstr(line, "Block");
         frame->name = (char *)malloc(strlen(strptr));
         strcpy(frame->name, strptr);
         frame->name[strlen(strptr)-1] = '\0'; // get rid of newline
         frame->index = ii;
         sscanf(strptr + 6, "%d", &frame->block_id);
         ii++;
         do_ties = 0;
         if(options->debug >= 15) 
             printf("  %s \n", frame->name);
         }

      if(strstr(line, "Conversion parameters")) {
         strptr = strstr(line, ":");
         strptr++;
         sscanf(strptr, "%lf %lf", &frame->cnvt_scale, &frame->min_pwr);
         }

   /* Get radiometric balancing coefficients */

      if(strstr(line, "Radiometric balancing offset")) {
           if(options->debug >= 20)
              printf("reading frame radiometric offset coefficients\n");
           frame->frm_offset = *get_coeffs(line, options);
           }

      if(strstr(line, "Radiometric balancing scale")) {
           if(options->debug >= 20)
              printf("reading frame radiometric scale coefficients\n");
           frame->frm_scale = *get_coeffs(line, options);
           }

    /* edge tie grunt */

      if(strstr(line, "Edge tie number")) {
         strptr = strstr(line, ":");
         strptr++;
         sscanf(strptr, "%d", &frame->frm_edgeties.n_ties);
         if(options->debug >= 15) 
            printf("number of frame ties %d\n", frame->frm_edgeties.n_ties);
         }

      if(strstr(line, "Edge tie spacing")) {
         strptr = strstr(line, ":");
         strptr++;
         sscanf(strptr, "%lf", &frame->frm_edgeties.spacing);
         if(options->debug >= 15) 
            printf("frame edgeties spacing %lf\n", frame->frm_edgeties.spacing);
         }

      if(strstr(line, "Edge ties") && frame->frm_edgeties.n_ties > 0) {
         frame->frm_edgeties.tie = (EdgeTie_t *) calloc(frame->frm_edgeties.n_ties, sizeof(EdgeTie_t));
         jj = 0;
         do_ties = 1;
         continue;
         }

      if(do_ties) {
         if(strlen(line) < 3) continue;
         tie = &frame->frm_edgeties.tie[jj];
         sscanf(line, "%lf %lf %lf", &tie->map_xy.x, &tie->map_xy.y, &tie->target);
         if(options->debug >= 20) 
            printf("%lf %lf %lf\n", tie->map_xy.x, tie->map_xy.y, tie->target);
         jj++;
         }
      }

   fclose(fp);

   /* Read Blocks.key file */

   sprintf(path, "%s/BLOCKS.KEY", INDEX_DIR);
   fp = fopen(path, "r");
   if(fp == NULL) {
      printf("missing %s\n", path);
      return 1;
      }

 /* count up the number of index frame references and allocate array of frames */

   ii = 0;
   while(fgets(line, 511, fp)) {
      if(strstr(line, "Block Index")) ii++;
      }
   rewind(fp);
   data->blocks = (Block_t *)calloc(ii, sizeof(Block_t ));
   data->n_blocks = ii;
   ii = 0;

   do_ties = 0;
   while(fgets(line, 511, fp)) {
      if(strstr(line, "Block Index")) {
         block = &data->blocks[ii];
         strptr = strstr(line, ":");
         strptr = strstr(strptr, "Block");
         sscanf(strptr, "%s", name);
         block->name = (char *)malloc(strlen(name));
         strcpy(block->name, name);
         sscanf(strptr + 6, "%d", &block->id);
         if(options->debug >= 15) 
             printf("  %s\n", block->name);
         ii++;
         do_ties = 0;
         }


   /* Get coefficients */

      if(strstr(line, "Radiometric balancing offset")) {
           if(options->debug >= 20)
              printf("reading block radiometric offset coefficients\n");
           block->blk_offset = *get_coeffs(line, options);
           }

      if(strstr(line, "Radiometric balancing scale")) {
           if(options->debug >= 20)
              printf("reading block radiometric scale coefficients\n");
           block->blk_scale = *get_coeffs(line, options);
           }

      if(strstr(line, "Geometric balancing parameters")) {
           if(options->debug >= 20)
              printf("reading block geometric coefficients\n");
           block->blk_geom = *get_coeffs(line, options);
           }

      /* ---- Read edge ties for block ---- */

      if(strstr(line, "Edge tie number")) {
         strptr = strstr(line, ":");
         strptr++;
         sscanf(strptr, "%d", &block->blk_edgeties.n_ties);
         if(options->debug >= 15) 
            printf("number of block ties %d\n", block->blk_edgeties.n_ties);
         }

      if(strstr(line, "Edge tie spacing")) {
         strptr = strstr(line, ":");
         strptr++;
         sscanf(strptr, "%lf", &block->blk_edgeties.spacing);
         if(options->debug >= 15) 
            printf("block edgeties spacing %lf\n", block->blk_edgeties.spacing);
         }

      if(strstr(line, "Edge ties") && block->blk_edgeties.n_ties > 0) {
         block->blk_edgeties.tie = (EdgeTie_t *) calloc(block->blk_edgeties.n_ties, sizeof(EdgeTie_t));
         jj = 0;
         do_ties = 1;
         continue;
         }

      if(do_ties) {
         if(strlen(line) < 3) continue;
         tie = &block->blk_edgeties.tie[jj];
         sscanf(line, "%lf %lf %lf", &tie->map_xy.x, &tie->map_xy.x, &tie->target);
         if(options->debug >= 20) 
            printf("%lf %lf %lf\n", tie->map_xy.x, tie->map_xy.x, tie->target);
         jj++;
         }
      }

   fclose(fp);

   calculate_output_parameters(data, options);

   return 0;
}

/*fs----------------------------------------------------------------------------

    Procedure:   get_coeffs

    Purpose:   extract a list of coefficients from a line

    Exits:   coefficients upon success, else null

----------------------------------------------------------------------------fe*/

coeffs_t *get_coeffs(char *line, Options_t *options)
{
   coeff_t *coeffptr = NULL;
   coeffs_t *coeffs;
   char *strptr = strstr(line, ":");
   char temp_line[512], *next_word;

   coeffs = (coeffs_t *)calloc(1, sizeof(coeffs_t));

   strptr++;
   strcpy(temp_line, strptr);

   next_word = strtok(temp_line, " ");

   while(next_word) {
      if(coeffs->firstcoeff == NULL) {
          coeffptr = coeffs->firstcoeff = (coeff_t *)calloc(1, sizeof(coeff_t));
      } else  {
          coeffptr->next = (coeff_t *)calloc(1, sizeof(coeff_t));
          coeffptr = coeffptr->next;
          }
      sscanf(next_word, "%lf", &coeffptr->value);
      if(options->debug >= 20)
          printf("     %le\n", coeffptr->value);
      coeffs->n_coeffs++;
      next_word = strtok(NULL, " ");
      if(next_word[0] =='\n') break;
      }
   return coeffs;
}

/*fs----------------------------------------------------------------------------

    Procedure:   calculate_output_parameters

    Purpose:   Calculate output file size and extents as well as
               the subtile pixel coordinates for each subtile

    Exits:   Exit status is 0 on success, 1 on failure

----------------------------------------------------------------------------fe*/

int calculate_output_parameters(Data_t *data, Options_t *options)
{
   Subtile_t *sub;
   Image_t *out, *index = NULL;
   int    ii, jj;
   double min_x, min_y, max_x, max_y;

   min_x = data->subs[0].min_x;
   min_y = data->subs[0].min_y;
   max_x = data->subs[0].max_x;
   max_y = data->subs[0].max_y;

 /* get overall extent */

   for(ii = 0; ii < data->n_subs; ii++) {
      sub = &data->subs[ii];
      if(min_x > sub->min_x) min_x = sub->min_x;
      if(min_y > sub->min_y) min_y = sub->min_y;
      if(max_x < sub->max_x) max_x = sub->max_x;
      if(max_y < sub->max_y) max_y = sub->max_y;
      }

   out = (Image_t *)calloc(1, sizeof(Image_t));

   if(options->index_file)
      index = (Image_t *)calloc(1, sizeof(Image_t));

   out->min_x = min_x;
   out->min_y = min_y;
   out->max_x = max_x;
   out->max_y = max_y;

   out->size_x = (max_x - min_x) / data->image_res;
   out->size_y = (max_y - min_y) / data->image_res;

   if(index) {
      index->size_x = (max_x - min_x) / data->index_res;
      index->size_y = (max_y - min_y) / data->index_res;
      }

   if(options->debug >= 5) {
        printf("output extent %.0lf to %.0lf Easting,  %.0lf to %.0lf Northing\n",
                min_x, max_x, min_y, max_y);
        printf("image size will be %d  %d\n", out->size_x, out->size_y);
        }

   for(ii = 0; ii < data->n_subs; ii++) {
      sub = &data->subs[ii];
      sub->img_ul_x = (sub->min_x - out->min_x) / data->image_res;
      sub->img_lr_x = (sub->max_x - out->min_x) / data->image_res;
      sub->img_ul_y = (out->max_y - sub->max_y) / data->image_res;
      sub->img_lr_y = (out->max_y - sub->min_y) / data->image_res;
      sub->index_ul_x = (sub->min_x - out->min_x) / data->index_res;
      sub->index_ul_y = (out->max_y - sub->max_y) / data->index_res;
      if(options->debug >= 20) 
          printf("%s:  %d  %d  %d  %d\n", sub->name, sub->img_ul_x, 
                  sub->img_lr_x, sub->img_ul_y, sub->img_lr_y);
      }

   for(ii = 0; ii < data->n_frames; ii++) {
      for(jj = 0; jj < data->n_blocks; jj++) {
         if(data->frames[ii].block_id == data->blocks[jj].id)
            break;
         }
      data->frames[ii].block = &data->blocks[jj];
      if(options->debug >= 25)
          printf("assigning block %d to frame %s\n", data->blocks[jj].id, data->frames[ii].name);
      }

   data->output_image = out;
   data->index_image =  index;

   return 0;
}

/*fs----------------------------------------------------------------------------

    Procedure:   calculate_sub

    Purpose:   Get sigma nought value for a subtile and write it to output

    Exits:   Exit status is 0 on success, 1 on failure

----------------------------------------------------------------------------fe*/

int calculate_sub(Subtile_t *sub, Options_t *options, Data_t *data)
{
   char *name = sub->name;
   FILE *fp;
   char path[256];
   static char fn[] = "calculate_sub";
   int   n_pixels = data->image_size;
   int buf_size = n_pixels * n_pixels;
   int scale = data->index_res / data->image_res;
   int i_size = data->index_size * data->index_size;

   int ii, jj, offset, i_offset;
   short *buf = NULL;
   unsigned char *i_buf = NULL; 

   float min_x, max_y;
   short *out_buf     = NULL;
   short no_data_val  = NO_DATA_VAL;
   short out_null     = OUT_NULL;
   float data_scale   = DATA_SCALE;
   short out_val, off = OFFSET;

   if(options->do_point) {
       if((options->map_x < sub->min_x) ||
         (options->map_x > sub->max_x) ||
         (options->map_y < sub->min_y) ||
         (options->map_y > sub->max_y)) return 0;
      }
 
   buf = (short *)calloc(buf_size, sizeof(short));
   out_buf = (short *)calloc(buf_size, sizeof(short));
   i_buf = (unsigned char *)calloc(i_size, 1);

   if(options->debug >= 1)
       printf("processing %s\n", sub->name);

   /* ---- Open <tile name>.IMG file ---- */
   sprintf(path,"IMAGES.DIR/%s.IMG", name);
   if ((fp = fopen(path, "rb")) == NULL) {
      printf("%s: Unable to open %s\n", fn, path);
      return 1;
      }

   fread(buf, sizeof(short), buf_size, fp);
   fclose(fp);

   /* ---- Open <tile name>.IDX file ---- */
   sprintf(path,"INDICES.DIR/%s.IDX", name);
   if ((fp = fopen(path, "rb")) == NULL) {
      printf("%s: Unable to open %s\n", fn, path);
      return 1;
     }

   fread(i_buf, 1, i_size, fp);

   /* ---- Close up ---- */
   fclose(fp);
   fp = NULL;

   min_x = sub->min_x;
   max_y = sub->max_y;

  /* do the point thing */
   if(options->do_point) {
      if(options->debug > 0) options->debug += 30; // crank up the debug
      ii = (sub->max_y - options->map_y) / data->image_res;
      jj = (options->map_x - sub->min_x) / data->image_res;
      offset = ii * n_pixels + jj;
      i_offset = ii/scale * n_pixels/scale+ jj/scale;
      data->value = (double)buf[offset];
      if (data->value == no_data_val) {
         printf("%lf %lf: %lf\n", options->map_x, options->map_y, (double)no_data_val);
         return 0;
         }
      data->index_value = (int)i_buf[i_offset];
      data->x = options->map_x;
      data->y = options->map_y;
      GetSigma0(options, data);
      out_val = (short)((data->s0 + off) * data_scale) - OUT_OFFSET;
      printf("%lf %lf: %lf\n", options->map_x, options->map_y, data->s0);
      if(options->debug > 40) // print this if the user specified a debug level >= 10
         printf("16bit encoded value = %d\n", out_val);
      return 0;
      }

   for(ii = 0; ii < n_pixels; ii++) {
      for(jj = 0; jj < n_pixels ; jj++) {
         offset = ii * n_pixels + jj;
         i_offset = ii/scale * n_pixels/scale+ jj/scale;
         data->value = (double)buf[offset];

      /* ---- Skip reading index image if no data value ---- */
         if (data->value == no_data_val) {
             out_buf[offset] = out_null;
             continue;
             }

       /* ---- Get value of pixel ---- */
       data->index_value = (int)i_buf[i_offset];
       if(data->index_value >= data->n_frames) {
           printf("%s: %d at %d %d exceeds index range %d\n", name, 
                 data->index_value, jj/scale, ii/scale, data->n_frames);
           exit(1);
           }

        data->x = min_x + jj * data->image_res;
        data->y = max_y - ii * data->image_res;

   /* ---- Convert value ---- */

        GetSigma0(options, data);
        if(data->s0 == no_data_val) {
            out_buf[offset] = out_null;
            continue;
            }
        if(data->s0 < -30) data->s0 = -30;
        if(data->s0 > 10) data->s0 = 10;
        out_buf[offset] = (short)((data->s0 + off) * data_scale) - OUT_OFFSET;
        }
   }

  sub->buf = out_buf; 
  sub->i_buf = i_buf; 
  if(write_sub(sub, data, options)) {
     printf("error writing subtile\n");
     return 1;
     }
  free(out_buf);
  free(i_buf);
  return 0;
}



/*fs----------------------------------------------------------------------------

    Procedure:   int GetSigma0 (options, data)

    Purpose:   To get value from line & sample that a selected x,y map to.

    Arguments:   Options_t *options - Pointer to command line options.
      Data_t   *data      - Pointer to program data structure.

  Other than getting the transformation parameters from the data structures,
  this was taken directly from Pete's code.
 
    Returns:   Returns 0 on success. Otherwise, returns 1 on failure.

----------------------------------------------------------------------------fe*/
int GetSigma0 (
          Options_t *options,
          Data_t *data)
{
  static char fn[] = "GetSigma0";
  double scale, offset, min, max;
  double x1, y1, aa, bb, cc, dd, diff;
  coeff_t *coeffptr;
  Frame_t *frame;
  Block_t *block;

  frame = &data->frames[data->index_value];
  block = frame->block;
  /* ---- Set the min max values alowed for the data type ---- */
  min = 0;
  max = 32767;

  /* ---- Initialize reversal value ---- */
  data->s0 = data->value;

  /* ---- Reverse edge balancing for block ---- */
  if (block->blk_edgeties.n_ties > 0) {
    offset = FindOffsetAtPt
      (data->x, data->y, block->blk_edgeties.tie,
       block->blk_edgeties.n_ties, block->blk_edgeties.spacing);
    data->s0 -= offset;
    if (data->s0 < min)
      data->s0 = min;
    else if (data->s0 > max)
      data->s0 = max;
    if (options->debug >= 30) {
      printf("%s: Reversing edge offset (%lf): %lf\n",
         fn, offset, data->s0);
    }
  }


  /* ---- Reverse grand_rad ---- */
  if (block->blk_offset.n_coeffs > 0) {
    if ((offset = ApplyEqnAtPt(data->x, data->y, &block->blk_offset)) == -9999) {
      printf("%s: Invalid equations\n", fn);
      return(1);
    }

    data->s0 -= offset;
    if (options->debug >= 30) {
      printf("%s: Reversing offset (%lf): %lf\n",
         fn, offset, data->s0);
    }

    if ((scale = ApplyEqnAtPt(data->x, data->y, &block->blk_scale)) == -9999) {
      printf("%s: Invalid equations\n", fn);
      return(1);
    }
    scale = pow(10, scale);
    data->s0 /= scale;
    if (options->debug >= 30) {
      printf("%s: Reversing scale (%lf): %lf\n",
         fn, scale, data->s0);
    }
  }

  /* >>>> Reverse grand_geo <<<< */

  /* ---- Init ---- */
  x1 = data->x;
  y1 = data->y;
  if (block->blk_geom.n_coeffs != 4) {
    printf("%s: Invalid geometric equation coefficients", fn);
    return(1);
  }
  coeffptr = block->blk_geom.firstcoeff;
  aa  = coeffptr->value;
  coeffptr = coeffptr->next;
  bb  = coeffptr->value;
  coeffptr = coeffptr->next;
  cc  = coeffptr->value;
  coeffptr = coeffptr->next;
  if (cc == 0) {
    printf("%s: Invalid geometric equation coefficients", fn);
    return(1);
  }
  dd  = coeffptr->value;
  diff = dd / cc;

  /* ---- Compute x, y location prior to geometric transformation ---- */
  data->geo_y = (diff * x1 + y1 - diff * aa + bb) / (diff * dd + cc);
  data->geo_x = (x1 - aa - dd * y1) / cc;

    if (options->debug >= 30) {
    printf("%s: Reversing geometric adj:", fn);
    printf(" (%.2lf, %.2lf)->(%.2lf, %.2lf)\n",
       data->x, data->y, data->geo_x, data->geo_y);
  }

  /* ---- Reverse edge balancing using coords from grand_geo reversal ---- */
  if (frame->frm_edgeties.n_ties > 0) {
    offset = FindOffsetAtPt(data->geo_x, data->geo_y, frame->frm_edgeties.tie,
            frame->frm_edgeties.n_ties, frame->frm_edgeties.spacing);
    data->s0 -= offset;
    if (options->debug >= 30) {
      printf("%s: Reversing edge offset (%lf): %lf\n",
         fn, offset, data->s0);
    }
  }

  /* ---- Reverse rad_bal using coords from grand_geo reversal ---- */
  if (frame->frm_offset.n_coeffs > 0) {
    if ((offset = ApplyEqnAtPt(data->geo_x, data->geo_y,
               &frame->frm_offset)) == -9999) {
      printf("%s: Invalid equation\n", fn);
      return(1);
    }
    data->s0 -= offset;
    if (options->debug >= 30) {
      printf("%s: Reversing offset (%lf): %lf\n",
         fn, offset, data->s0);
    }

    if ((scale = ApplyEqnAtPt(data->geo_x, data->geo_y,
               &frame->frm_scale)) == -9999) {
      printf("%s: Invalid equation", fn);
      return(1);
    }
    scale = pow(10, scale);
    data->s0 /= scale;
    if (options->debug >= 30) {
      printf("%s: Reversing scale (%lf): %lf\n", fn, scale, data->s0);
    }
  }
/* this is the mamm version of the scale/offset calculation */
#ifdef MAMM

  data->s0 -= frame->min_pwr;

    if (options->debug >= 30) {
    printf("%s: Reversing offset (%lf): %lf\n",
            fn, frame->min_pwr, data->s0);
  }

  data->s0 /= frame->cnvt_scale;

    if (options->debug >= 30) {
    printf("%s: Reversing scale (%lf): %lf\n",
       fn, frame->cnvt_scale, data->s0);
  }

#else

/* this is the amm1 version.  */

  data->s0 /= frame->cnvt_scale;

  /* ---- Reverse original scale / offset ---- */
    if (options->debug >= 30) {
    printf("%s: Reversing scale (%lf): %lf\n",
       fn, frame->cnvt_scale, data->s0);
  }

  data->s0 -= frame->min_pwr;

    if (options->debug >= 30) {
    printf("%s: Reversing offset (%lf): %lf\n",
            fn, frame->min_pwr, data->s0);
  }

#endif

  /* ---- Reverse conversion to amplitude from power ---- */
    if (options->debug >= 30) {
    printf("%s: Amplitude: %lf\n", fn, data->s0);
  }
  data->s0 *= data->s0;
    if (options->debug >= 30) {
    printf("%s: Power: %lf\n", fn, data->s0);
  }

  /* ---- Obtain sigma0 from power ---- */
  data->s0 = 10 * log10(data->s0);
    if (options->debug >= 30) {
    printf("%s: Sigma Nought: %lf\n", fn, data->s0);
  }

  /* ---- Return success ---- */
  return(0);
}

/*fs----------------------------------------------------------------------------

   Procedure: double FindOffsetAtPt(xx, yy, tie_array, n_ties, spacing)

   Purpose:   To compute offset at point based on distance to tie points and
      corresponding tie offset values.

   Arguments:
      double       xx, yy     - Map coordinates at which to compute offset
      EdgeTie_t    *tie_array - Ptr to array of tie points
      int          n_ties     - Number of ties in tie_array
      double       spacing    - Tie point spacing

   Returns:   Offset value.

----------------------------------------------------------------------------fe*/
double FindOffsetAtPt(
   double     xx,
   double     yy,
   EdgeTie_t  *tie_array,
   int        n_ties,
   double     spacing)
{
   int max, min, mid, ii;
   double mid_x;
   double dx, dy, dist;
   double offset;


   /* ---- Init ---- */
   offset = 0;
   min = 0;
   max = n_ties - 1;


   /* >>>> Find ties within spacing distance <<<< */

   /* ---- Find lowest valid x value ---- */
   while (tie_array[max].map_xy.x > tie_array[min].map_xy.x)
   {
      mid = (max + min) / 2;
      mid_x = tie_array[mid].map_xy.x;
      if (mid_x >= xx - spacing) {
         max = mid;
      }
      else {
     if (min == mid) {
            min = max;
     }
         else {
            min = mid;
     }
      }
   }

   /* ---- Step through ties until out of target range ---- */
   for (ii = min; ii < n_ties; ii++) {

      /* ---- See if we've passed x range ---- */
      dx = tie_array[ii].map_xy.x - xx;
      if (dx > spacing) {
         break;
      }

      /* ---- See if y is within range ---- */
      dy = tie_array[ii].map_xy.y - yy;
      if (dy < 0) dy = -dy;
      if (dy > spacing) continue;

      /* ---- Check distance ---- */
      dist = sqrt(dx * dx + dy * dy);
      if (dist > spacing) continue;

      /* ---- See if point is close enough to contribute to offset ---- */
      if (dist < spacing) {

         /* ---- Pt found, calculate contribution to offset ---- */
         offset += (tie_array[ii].target - tie_array[ii].avg) *
               ((spacing - dist) / spacing);

      }
   }

   /* ---- Return offset ---- */
   return(offset);
}

/*fs----------------------------------------------------------------------------

   Procedure:   double ApplyEqnAtPt(xx, yy, coeffs)

   Purpose: To apply equation to a map point.

   Arguments:
      double xx, yy   - Map coords of location to apply equation.
      coeff_t *coeffs - Structure containing equation coefficients.

   Returns: Returns value on success or -9999 on failure.

----------------------------------------------------------------------------fe*/
double ApplyEqnAtPt(
   double xx,
   double yy,
   coeffs_t *coeffs)
{
   static char fn[] = "ApplyEqnAtPt";
   int order;
   int terms;
   int ii;
   double pow_x, pow_y, offset;
   coeff_t *coeffptr;

   /* ---- Init ---- */
   order  = 0;
   terms  = 0;
   offset = 0;
   if ((coeffs == NULL) ||
      ((coeffptr = coeffs->firstcoeff) == NULL)) {
      printf("%s: Invalid list of equation coefficients\n", fn);
      return -9999;
   }

   /* ---- Get order from number of terms ---- */
   while (terms < coeffs->n_coeffs ) {
      order++;
      terms = 0;
      for (ii = order; ii > 0; ii-- ) {
         terms += ii;
      }
      terms += order + 1;
   }

   /* ---- Calculate value for x, y ---- */
   while (order > 0) {
      pow_y = 0;
      for (pow_x = order; pow_x >= 0; pow_x--) {
         offset += pow(xx, pow_x) * pow(yy, pow_y) * coeffptr->value;
         pow_y++;
         coeffptr = coeffptr->next;

         /* ---- Check list just in case ---- */
         if (coeffptr == NULL) {
            printf("%s: Invalid list of equation coefficients\n", fn);
            return -9999;
     }
      }
      order--;
   }
   offset += coeffptr->value;

   /* ---- Return value ---- */
   return offset;
}

/*fs----------------------------------------------------------------------------

   Procedure:   write_sub

   Purpose:     write the data from a subtile into the output 
                data and index files.

----------------------------------------------------------------------------fe*/

int write_sub(Subtile_t *sub, Data_t *data, Options_t *options)
{
   long int offset;
   int jj, ii, out_i;
   void *data_ptr;
   int  n_bytes = data->image_size * sizeof(short);
   if(options->debug > 0)
       printf("writing %s\n", sub->name);

   for(ii = 0; ii < data->image_size; ii++) {
      out_i = sub->img_ul_y + ii;
      offset = (out_i * data->output_image->size_x + sub->img_ul_x) * sizeof(short);
      lseek(data->output_image->fd, offset, SEEK_SET);
      jj = ii * data->image_size;
      data_ptr = (void *)&sub->buf[jj];
      write(data->output_image->fd, data_ptr, n_bytes); 
      }

   if(!data->index_image)
      return 0;

   n_bytes = data->index_size;
   for(ii = 0; ii < data->index_size; ii++) {
      out_i = sub->index_ul_y + ii;
      offset = (out_i * data->index_image->size_x + sub->index_ul_x);
      lseek(data->index_image->fd, offset, SEEK_SET);
      jj = ii * data->index_size;
      data_ptr = (void *)&sub->i_buf[jj];
      write(data->index_image->fd, data_ptr, n_bytes); 
      }

   return 0;
}

/*fs----------------------------------------------------------------------------

   Procedure:   give_head

   Purpose:     Generate corners files and rams format header files for
                both the data and index output.  

----------------------------------------------------------------------------fe*/

int give_head(Data_t *data, Options_t *options)
{
   char path[1024];
   FILE *fp;

   sprintf(path, "%s.h", options->output_file);
   fp = fopen(path, "w");

   fprintf(fp, "#       @(#)image.c     2.2  7/13/0\n");
   fprintf(fp, "#       (c) Copyright 2005 Vexcel Corporation.\n");
   fprintf(fp, "#       All Rights Reserved\n");
   fprintf(fp, "#\n");
   fprintf(fp, "Vexcel Image Header\n");
   fprintf(fp, "lines   %d\n", data->output_image->size_y);
   fprintf(fp, "pixels  %d\n", data->output_image->size_x);
   fprintf(fp, "banding BIL\n");
   fprintf(fp, "bands   1\n");
   fprintf(fp, "data    short\n");
   fprintf(fp, "endian  BIG\n");
   fprintf(fp, "file    '%s'\n", options->output_file);
   fclose(fp);

   sprintf(path, "%s.corners", options->output_file);
   fp = fopen(path, "w");

   fprintf(fp, "#       @(#)rectxy.c    2.2     7/13/0\n");
   fprintf(fp, "#       (c) Copyright 2005 Vexcel Corporation.\n");
   fprintf(fp, "#       All Rights Reserved\n");
   fprintf(fp, "#\n");
   fprintf(fp, "Vexcel Rectangular XY coordinates\n");
   fprintf(fp, "max_x   %f\n", (float)data->output_image->max_x);
   fprintf(fp, "max_y   %f\n", (float)data->output_image->max_y);
   fprintf(fp, "min_x   %f\n", (float)data->output_image->min_x);
   fprintf(fp, "min_y   %f\n", (float)data->output_image->min_y);

   fclose(fp);

   if(options->index_file == NULL)
      return 0;

   sprintf(path, "%s.h", options->index_file);
   fp = fopen(path, "w");

   fprintf(fp, "#       @(#)image.c     2.2  7/13/0\n");
   fprintf(fp, "#       (c) Copyright 2005 Vexcel Corporation.\n");
   fprintf(fp, "#       All Rights Reserved\n");
   fprintf(fp, "#\n");
   fprintf(fp, "Vexcel Image Header\n");
   fprintf(fp, "lines   %d\n", data->index_image->size_y);
   fprintf(fp, "pixels  %d\n", data->index_image->size_x);
   fprintf(fp, "banding BIL\n");
   fprintf(fp, "bands   1\n");
   fprintf(fp, "data    byte\n");
   fprintf(fp, "endian  BIG\n");
   fprintf(fp, "file    '%s'\n", options->index_file);
   fclose(fp);

   return 0;
}

/*fs----------------------------------------------------------------------------

   Procedure:   prepare_output

   Purpose:     Generate output files for both the index and 
                data initialized with no_data values.  
                the file descriptors for each stay open and
                are returned with the images. 

----------------------------------------------------------------------------fe*/

int prepare_output(Data_t *data, Options_t *options)
{
   int fd;
   int ii, n_bytes;
   void *buf;
   short *s_buf;
   unsigned char *i_buf;

   if(options->debug >= 1)
       printf("preparing output image\n");

   n_bytes = data->output_image->size_x * sizeof(short);
   buf = calloc(data->output_image->size_x, sizeof(short));
   s_buf = (short *)buf;
   for(ii = 0; ii <  data->output_image->size_x; ii++) 
      s_buf[ii] = -32767;

   fd = open(options->output_file, O_WRONLY | O_CREAT | O_TRUNC, 0664);

   for(ii = 0; ii < data->output_image->size_y; ii++) {
       write(fd, buf, n_bytes); 
       }
   data->output_image->fd = fd;
   free(buf);

   if(options->index_file == NULL)
      return 0;

   n_bytes = data->index_image->size_x;
   buf = calloc(data->index_image->size_x, 1);
   i_buf = (unsigned char *)buf;
   for(ii = 0; ii <  data->index_image->size_x; ii++) 
      i_buf[ii] = 255;
   fd = open(options->index_file, O_WRONLY | O_CREAT | O_TRUNC, 0664);

   for(ii = 0; ii < data->index_image->size_y; ii++) {
       write(fd, buf, n_bytes);
       }
   data->index_image->fd = fd;
   free(buf);

   return 0;
}
