/* Standard Library Includes */
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* IAS Library Includes */
#include "ias_angle_gen_includes.h"
#include "ias_logging.h"
#include "ias_odl.h"      
#include "ias_miscellaneous.h"  
#include "ias_angle_gen_distro.h"
#include "ias_angle_gen_private.h"

/* Local Defines */
#define MAX_FIELD_NAME  32 /* Maximum number of characters in an ODL field */
#define SECS_TOLERANCE .0000015 /* Epoch seconds tolerance */ 

/*******************************************************************************
 Name: get_odl_field
 
 Purpose: Retrieves the ODL field and does error checking on that field.
 
 Returns:
    Type = integer 
    SUCCESS / ERROR
 ******************************************************************************/
static int get_odl_field 
(
    int memory_size,       /* I: Size of memory to retrieve */
    IAS_ODL_TYPE type,     /* I: Type of ODL field to retreive */
    IAS_OBJ_DESC *tree,    /* I: ODL tree to retrieve the field */
    const char *group,     /* I: ODL group where field is located */
    const char *field,     /* I: ODL field name */
    int expected_count,    /* I: Expected values to retrieve */
    void *address          /* O: Pointer to storage address */
)
{
    int count;  /* Number of values retrieved */

    /* Retrieve the ODL field */
    if (ias_odl_get_field(address, memory_size, type, tree, group, field, 
        &count) != SUCCESS)
    {
        IAS_LOG_ERROR("Reading %s from the ODL tree", field);
        return ERROR;
    }
    
    /* Compare the actual and expected count */
    if (count != expected_count)
    {
        IAS_LOG_ERROR("Retrieved %d values expected values is %d for %s", 
            count, expected_count, field);
        return ERROR;
    }

    return SUCCESS;
}

/*******************************************************************************
 Name: construct_and_get_odl_field
 
 Purpose: Constructs the ODL field name and then retrieves the odl field.
 
 Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
static int construct_and_get_odl_field
(
    int memory_size,           /* I: Size of memory to retrieve */
    IAS_ODL_TYPE type,         /* I: Type of ODL field to retreive */
    IAS_OBJ_DESC *tree,        /* I: ODL tree to retrieve the field */
    const char *group,         /* I: ODL group where field is located */
    const char *field,         /* I: ODL field name */
    const char *field_prepend, /* I: Value prepended to field */
    int expected_count,        /* I: Expected values to retrieve */
    void *address              /* O: Pointer to storage address */
)
{
    char field_name[MAX_FIELD_NAME]; /* Field name */
    int status; 

    /* Read the number of L1T samples for this band */
    status = snprintf(field_name, sizeof(field_name), "%s_%s", field_prepend, 
        field);
    if (status < 0 || status >= sizeof(field_name))
    {
        IAS_LOG_ERROR("Constructing the odl field for %s_%s", field_prepend,
            field);
        return ERROR;
    }

    /* Retrieve the ODL field */
    if (get_odl_field(memory_size, type, tree, group, field_name, 
        expected_count, address) != SUCCESS)
    {
        return ERROR;
    }

    return SUCCESS;
}

/*******************************************************************************
Name: ias_angle_gen_read_ang_header

Purpose: Reads the header group from the ANG metadata file.

Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
static int ias_angle_gen_read_ang_header
(
    IAS_OBJ_DESC *odl_data,          /* I: Metadata ODL object */
    IAS_ANGLE_GEN_METADATA *metadata /* O: Metadata structure to load */
)       
{
    int index;                      /* Loop index */
    int band_list[IAS_MAX_NBANDS];  /* Local array for loading band  numbers */

    /* Read landsat scene id */
    if (get_odl_field(sizeof(metadata->landsat_scene_id), IAS_ODL_String, 
        odl_data, "FILE_HEADER", "LANDSAT_SCENE_ID", 1,
        &metadata->landsat_scene_id) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the satellite name */
    if (get_odl_field(sizeof(metadata->spacecraft_id), IAS_ODL_String, odl_data,
        "FILE_HEADER", "SPACECRAFT_ID", 1, &metadata->spacecraft_id) != SUCCESS)
    {
        return ERROR;
    }

    /* Read number of bands */
    if (get_odl_field(sizeof(metadata->num_bands), IAS_ODL_Int, odl_data, 
        "FILE_HEADER", "NUMBER_OF_BANDS", 1, &metadata->num_bands) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the band list */
    if (get_odl_field(metadata->num_bands * sizeof(metadata->num_bands), 
        IAS_ODL_Int, odl_data, "FILE_HEADER", "BAND_LIST", metadata->num_bands, 
        band_list) != SUCCESS)
    {
        return ERROR;
    }

    /* Store the band list numbers in the data structure */
    for (index = 0; index < metadata->num_bands; index++) 
    {
        int band_index = ias_sat_attr_convert_band_number_to_index(
            band_list[index]);
        if (band_index == ERROR)
        {
            IAS_LOG_ERROR("Retrieving the band index from the satellite "
                "attributes library for band number %d", band_list[index]);
            return ERROR;
        }
        
        /* Check that band index is in the valid range */
        if (band_index < 0 || band_index >= IAS_MAX_NBANDS)
        {
            IAS_LOG_ERROR("Band index %d is not in the valid index range",
                band_index);
            return ERROR;
        }

        metadata->band_metadata[band_index].band_number = band_list[index];
        metadata->band_present[band_index] = TRUE;
    }

    return SUCCESS;
}

/*******************************************************************************
Name: ias_angle_gen_read_ang_projection

Purpose: Reads the projection group from the ANG metadata file.

Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
static int ias_angle_gen_read_ang_projection
(
    IAS_OBJ_DESC *odl_data,          /* I: Metadata ODL object */
    IAS_ANGLE_GEN_METADATA *metadata /* O: Metadata structure to load */
)        
{
    double coords[2];      /* Local array for loading coordinate pairs  */
    char map_projection[4];/* Local map projection */    
    char ellipsoid[6];     /* Local ellipsoid */

    /* Read the ellipsoid axes*/
    if (get_odl_field(sizeof(coords), IAS_ODL_Double, odl_data, "PROJECTION",
        "ELLIPSOID_AXES", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->wgs84_major_axis = coords[0];
    metadata->wgs84_minor_axis = coords[1];

    /* Read the map projection */
    if (get_odl_field(sizeof(map_projection), IAS_ODL_String, odl_data, 
        "PROJECTION", "MAP_PROJECTION", 1, map_projection) != SUCCESS)
    {
        return ERROR;
    }

    /* Process the map projection */
    ias_misc_convert_to_uppercase(map_projection);
    if (strcmp(map_projection, "UTM") == 0)
    {
        metadata->projection.proj_code = UTM;
    }
    else if (strcmp(map_projection, "PS") == 0)
    {
        metadata->projection.proj_code = PS;
    }
    else if (strcmp(map_projection, "AEA") == 0)
    {
        metadata->projection.proj_code = ALBERS;
    }
    else
    {
        IAS_LOG_ERROR("MAP_PROJECTION parameter is invalid must be "
            "AEA, UTM or PS");
        return ERROR;
    }

    /* Read the projection units */
    if (get_odl_field(sizeof(metadata->units), IAS_ODL_String, odl_data, 
        "PROJECTION", "PROJECTION_UNITS", 1, metadata->units) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the projection datum */
    if (get_odl_field(sizeof(metadata->datum), IAS_ODL_String, odl_data, 
        "PROJECTION", "DATUM", 1, metadata->datum) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the projection spheroid */
    if (get_odl_field(sizeof(ellipsoid), IAS_ODL_String, odl_data, "PROJECTION",
        "ELLIPSOID", 1, ellipsoid) != SUCCESS)
    {
        return ERROR;
    }
    
    metadata->projection.spheroid = WGS84_SPHEROID;
    ias_misc_convert_to_uppercase(ellipsoid);
    if (strcmp(ellipsoid, "WGS84") != 0)
    {
        IAS_LOG_ERROR("Unexpected ELLIPSOID: %s", ellipsoid);
        return ERROR;
    }

    /* Read the projection zone only if projection is UTM */
    if (metadata->projection.proj_code == UTM)
    {
        if (get_odl_field(sizeof(&metadata->projection.zone), IAS_ODL_Int, 
            odl_data, "PROJECTION", "UTM_ZONE", 1, &metadata->projection.zone) 
            != SUCCESS)
        {
            return ERROR;
        }

        /* Validate the UTM zone */
        if (metadata->projection.zone < 1 || metadata->projection.zone > 60)
        {
            IAS_LOG_ERROR("UTM_ZONE must be (1 - 60)");
            return ERROR;
        }
    }

    /* Read in all the projection parameters */
    if (get_odl_field(sizeof(metadata->projection.parameters), IAS_ODL_Double,
        odl_data, "PROJECTION", "PROJECTION_PARAMETERS", IAS_PROJ_PARAM_SIZE,
        metadata->projection.parameters) != SUCCESS)
    {
        return ERROR;
    }

    /* Read in the upper left corner */
    if (get_odl_field(sizeof(coords), IAS_ODL_Double, odl_data, "PROJECTION", 
        "UL_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.upleft.x = coords[0];
    metadata->corners.upleft.y = coords[1];

    /* Read in the upper right corner */
    if (get_odl_field(sizeof(coords), IAS_ODL_Double, odl_data, "PROJECTION", 
        "UR_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.upright.x = coords[0];
    metadata->corners.upright.y = coords[1];

    /* Read in the lower left corner */
    if (get_odl_field(sizeof(coords), IAS_ODL_Double, odl_data,
        "PROJECTION", "LL_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.loleft.x = coords[0];
    metadata->corners.loleft.y = coords[1];

    /* Read in the lower right corner */
    if (get_odl_field(sizeof(coords), IAS_ODL_Double, odl_data,
        "PROJECTION", "LR_CORNER", 2, coords) != SUCCESS)
    {
        return ERROR;
    }
    metadata->corners.loright.x = coords[0];
    metadata->corners.loright.y = coords[1];

    return SUCCESS;
}

/*******************************************************************************
Name: ias_angle_gen_read_ang_ephemeris

Purpose: Reads the ephemeris group from the ANG metadata file.

Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
static int ias_angle_gen_read_ang_ephemeris
(
    IAS_OBJ_DESC *odl_data,          /* I: Metadata ODL object */
    IAS_ANGLE_GEN_METADATA *metadata /* O: Metadata structure to load */
)        
{
    int index;          /* Loop index */
    int year;           /* Year of acquisition */
    int day;            /* Day of acquisition */
    double seconds;     /* Seconds of acquisition */
    int solar_points;   /* Number of solar points */
    size_t dbuffer_size;/* Size of dbuffer array */
    double *dbuffer;    /* Local array for loading values */

    /* Read the epoch year */
    if (get_odl_field(sizeof(year), IAS_ODL_Int, odl_data, "EPHEMERIS",
        "EPHEMERIS_EPOCH_YEAR", 1, &year) != SUCCESS)
    {
        return ERROR;
    }
    metadata->ephem_epoch_time.year = year;
    
    /* Ephemeris epoch time - day of year */
    if (get_odl_field(sizeof(day), IAS_ODL_Int, odl_data, "EPHEMERIS",
        "EPHEMERIS_EPOCH_DAY", 1, &day) != SUCCESS)
    {
        return ERROR;
    }
    metadata->ephem_epoch_time.day = day;

    /* Ephemeris epoch time - second of day */
    if (get_odl_field(sizeof(metadata->ephem_epoch_time.seconds), 
        IAS_ODL_Double, odl_data, "EPHEMERIS", "EPHEMERIS_EPOCH_SECONDS", 1,
        &metadata->ephem_epoch_time.seconds) != SUCCESS)
    {
        return ERROR;
    }

    /* Number of ephemeris points */
    if (get_odl_field(sizeof(metadata->ephem_count), IAS_ODL_Int, odl_data, 
        "EPHEMERIS", "NUMBER_OF_POINTS", 1, &metadata->ephem_count) != SUCCESS)
    {
        return ERROR;
    }

    /* Allocate space for the input buffer */
    dbuffer_size = metadata->ephem_count * sizeof(double);
    dbuffer = (double *)malloc(dbuffer_size);
    if (!dbuffer)
    {
        IAS_LOG_ERROR("Allocating ephemeris buffer");
        return ERROR;
    }
    
    /* Initialize the metadata structure */
    if (ias_angle_gen_initialize(metadata) != SUCCESS)
    {
        IAS_LOG_ERROR("Allocating ephemeris structures");
        free(dbuffer);
        return ERROR;
    }

    /* Ephemeris sample time offsets from epoch */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_TIME", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    { 
        metadata->ephemeris[index].sample_time = dbuffer[index];
    }

    /* Ephemeris ECEF X coordinates */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_ECEF_X", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->ephemeris[index].position.x = dbuffer[index];
    }

    /* Ephemeris ECEF Y coordinates */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_ECEF_Y", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    {
        metadata->ephemeris[index].position.y = dbuffer[index];
    }

    /* Ephemeris ECEF Z coordinates */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "EPHEMERIS", 
        "EPHEMERIS_ECEF_Z", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->ephemeris[index].position.z = dbuffer[index];
    }

    /* Solar vector data group */
    /* Read the Solar Vector epoch year */
    if (get_odl_field(sizeof(year), IAS_ODL_Int, odl_data, "SOLAR_VECTOR", 
        "SOLAR_EPOCH_YEAR", 1, &year) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    if ((int)metadata->ephem_epoch_time.year != year)
    {
        IAS_LOG_ERROR("Solar epoch year %d and Ephemeris epoch year %d are not"
            "equal", year, (int)metadata->ephem_epoch_time.year);
        free(dbuffer);
        return ERROR;
    }
    
    /* Solar vector epoch time - day of year */
    if (get_odl_field(sizeof(day), IAS_ODL_Int, odl_data, "SOLAR_VECTOR",
        "SOLAR_EPOCH_DAY", 1, &day) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }
    
    if ((int)metadata->ephem_epoch_time.day != day)
    {
        IAS_LOG_ERROR("Solar epoch day %d and Ephemeris epoch day %d are not"
            "equal", day, (int)metadata->ephem_epoch_time.day);
        free(dbuffer);
        return ERROR;
    }

    /* Solar vector epoch time - seconds of day */
    if (get_odl_field(sizeof(seconds), IAS_ODL_Double, odl_data, "SOLAR_VECTOR",
        "SOLAR_EPOCH_SECONDS", 1, &seconds) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    if (fabs(metadata->ephem_epoch_time.seconds - seconds) > SECS_TOLERANCE)
    {
        IAS_LOG_ERROR("Solar epoch seconds %lf and Ephemeris epoch seconds "
            "%lf difference is higher than %lf tolerance", seconds, 
            metadata->ephem_epoch_time.seconds, SECS_TOLERANCE);
        free(dbuffer);
        return ERROR;
    }

    /* Earth to Sun distance */
    if (get_odl_field(sizeof(metadata->earth_sun_distance), IAS_ODL_Double, 
        odl_data, "SOLAR_VECTOR", "EARTH_SUN_DISTANCE", 1, 
        &metadata->earth_sun_distance) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    /* Number of solar vector points */
    if (get_odl_field(sizeof(solar_points), IAS_ODL_Int, odl_data, 
        "SOLAR_VECTOR", "NUMBER_OF_POINTS", 1, &solar_points) != SUCCESS)
    {
        free(dbuffer);
        return ERROR;
    }

    if (metadata->ephem_count != solar_points)
    {
        IAS_LOG_ERROR("Solar points %d and Ephemeris points %d are not"
            "equal", solar_points, metadata->ephem_count);
        free(dbuffer);
        return ERROR;
    }

    /* Solar vector sample time offsets from epoch */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "SOLAR_VECTOR", 
        "SAMPLE_TIME", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->solar_vector[index].sample_time = dbuffer[index];
    }

    /* Solar vector ECEF X coordinates */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "SOLAR_VECTOR", 
        "SOLAR_ECEF_X", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    { 
        metadata->solar_vector[index].position.x = dbuffer[index];
    }

    /* Solar vector ECEF Y coordinates */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "SOLAR_VECTOR",
        "SOLAR_ECEF_Y", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++)
    {
        metadata->solar_vector[index].position.y = dbuffer[index];
    }

    /* Solar vector ECEF Z coordinates */
    if (get_odl_field(dbuffer_size, IAS_ODL_Double, odl_data, "SOLAR_VECTOR",
        "SOLAR_ECEF_Z", metadata->ephem_count, dbuffer) != SUCCESS)
    {
        ias_angle_gen_free(metadata);
        free(dbuffer);
        return ERROR;
    }

    for (index = 0; index < metadata->ephem_count; index++) 
    {
        metadata->solar_vector[index].position.z = dbuffer[index];
    }

    /* Free the input buffer */
    free(dbuffer);

    return SUCCESS;
}

/*******************************************************************************
Name: ias_angle_gen_read_ang_sca

Purpose: Reads the SCA group from the ANG metadata file.

Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
static int ias_angle_gen_read_ang_sca
(
    IAS_OBJ_DESC *odl_data,               /* I: ODL object */
    int band_number,                      /* I: Band number */
    int sca_number,                       /* I: SCA number */
    char *band_group,                     /* I: Band ODL group to read from */
    IAS_ANGLE_GEN_IMAGE_RPC *sca_metadata /* O: SCA metadata to load */
)
{
    char band_sca_string[13];   /* Current Band SCA combination string */
    double offsets[2];          /* Local offsets */
    int status;                 /* Routine return placeholder */

    /* Assign the sca number */
    sca_metadata->sca_number = sca_number;

    /* Construct the band string */
    status = snprintf(band_sca_string, sizeof(band_sca_string), 
        "BAND%02d_SCA%02d", band_number, sca_number);
    if (status < 0 || status >= sizeof(band_sca_string))
    {
        IAS_LOG_ERROR("Constructing the band/sca string for SCA %d", 
            sca_number);
        return ERROR;
    }

    /* Mean height */
    if (construct_and_get_odl_field(sizeof(sca_metadata->mean_height), 
        IAS_ODL_Double, odl_data, band_group, "MEAN_HEIGHT", band_sca_string, 1,
        &sca_metadata->mean_height) != SUCCESS)
    {
        return ERROR;
    }

    /* Mean offset L1R line/samp */
    if (construct_and_get_odl_field(sizeof(offsets), IAS_ODL_Double, odl_data, 
        band_group, "MEAN_L1R_LINE_SAMP", band_sca_string, 2, offsets)
        != SUCCESS)
    {
        return ERROR;
    }
    sca_metadata->line_terms.l1r_mean_offset = offsets[0];        
    sca_metadata->samp_terms.l1r_mean_offset = offsets[1]; 

    /* Mean offset L1T line/samp */ 
    if (construct_and_get_odl_field(sizeof(offsets), IAS_ODL_Double, odl_data, 
        band_group, "MEAN_L1T_LINE_SAMP", band_sca_string, 2, offsets) 
        != SUCCESS)
    {
        return ERROR;
    }
    sca_metadata->line_terms.l1t_mean_offset = offsets[0];        
    sca_metadata->samp_terms.l1t_mean_offset = offsets[1];

    /* Line numerator coefficients */
    if (construct_and_get_odl_field(IAS_ANGLE_GEN_NUM_RPC_COEF * sizeof(double),
        IAS_ODL_Double, odl_data, band_group, "LINE_NUM_COEF", band_sca_string,
        IAS_ANGLE_GEN_NUM_RPC_COEF, sca_metadata->line_terms.numerator) 
        != SUCCESS)
    {
        return ERROR;
    }

    /* Line denominator coefficients */
    if (construct_and_get_odl_field((IAS_ANGLE_GEN_NUM_RPC_COEF - 1) 
        * sizeof(double), IAS_ODL_Double, odl_data, band_group, "LINE_DEN_COEF",
        band_sca_string, IAS_ANGLE_GEN_NUM_RPC_COEF - 1, 
        sca_metadata->line_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }
    
    /* Sample numerator coefficients */
    if (construct_and_get_odl_field(IAS_ANGLE_GEN_NUM_RPC_COEF * sizeof(double),
        IAS_ODL_Double, odl_data, band_group, "SAMP_NUM_COEF", band_sca_string,
        IAS_ANGLE_GEN_NUM_RPC_COEF, sca_metadata->samp_terms.numerator) 
        != SUCCESS)
    {
        return ERROR;
    }

    /* Sample denominator coefficients */
    if (construct_and_get_odl_field((IAS_ANGLE_GEN_NUM_RPC_COEF - 1) 
        * sizeof(double), IAS_ODL_Double, odl_data, band_group, "SAMP_DEN_COEF",
        band_sca_string, IAS_ANGLE_GEN_NUM_RPC_COEF - 1, 
        sca_metadata->samp_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    return SUCCESS;
}   

/*******************************************************************************
Name: ias_angle_gen_read_ang_band

Purpose: Reads the band group from the ANG metadata file.

Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
static int ias_angle_gen_read_ang_band
(
    IAS_OBJ_DESC *odl_data,           /* I: ODL object */
    IAS_ANGLE_GEN_BAND *band_metadata /* O: Metadata band structure to load */
)       
{
    int sca;                        /* Sca loop variable*/
    int status;                     /* Status placeholder */
    int sca_list[IAS_MAX_NSCAS];    /* Input buffer for SCA list */
    char band_group[11];            /* Band group name */
    double offsets[2];              /* Local offsets */
    char band_string[7];            /* Current band string */
    size_t num_size;                /* Numerator data size */
    size_t den_size;                /* Denominator data size */

    /* Initialize the numerator and denominator data sizes */
    num_size = IAS_ANGLE_GEN_ANG_RPC_COEF * sizeof(double);
    den_size = (IAS_ANGLE_GEN_ANG_RPC_COEF - 1) * sizeof(double);

    /* Construct the band string */
    status = snprintf(band_string, sizeof(band_string), "BAND%02d", 
        band_metadata->band_number);
    if (status < 0 || status >= sizeof(band_string))
    {
        IAS_LOG_ERROR("Constructing the band string");
        return ERROR;
    }

    /* Read data */
    /* Band information */
    /* Construct the current band's group name */
    status = snprintf(band_group, sizeof(band_group), "RPC_%s", band_string);
    if (status < 0 || status >= sizeof(band_group))
    {
        IAS_LOG_ERROR("Constructing the band group");
        return ERROR;
    }

    /* Read the number of SCAs for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->num_scas),
        IAS_ODL_Int, odl_data, band_group, "NUMBER_OF_SCAS", band_string, 1,
        &band_metadata->num_scas) != SUCCESS)
    {
        return ERROR;
    }

    if (construct_and_get_odl_field(sizeof(band_metadata->l1t_lines),
        IAS_ODL_Int, odl_data, band_group, "NUM_L1T_LINES", band_string, 1,
        &band_metadata->l1t_lines) != SUCCESS)
    {
        return ERROR;
    }
    
    /* Read the number of L1T samples for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->l1t_samps),
        IAS_ODL_Int, odl_data, band_group, "NUM_L1T_SAMPS", band_string, 1,
        &band_metadata->l1t_samps) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the number of L1R lines for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->l1r_lines),
        IAS_ODL_Int, odl_data, band_group, "NUM_L1R_LINES", band_string, 1,
        &band_metadata->l1r_lines) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the number of L1R samples for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->l1r_samps),
        IAS_ODL_Int, odl_data, band_group, "NUM_L1R_SAMPS", band_string, 1,
        &band_metadata->l1r_samps) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the pixel size for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->pixel_size),
        IAS_ODL_Double, odl_data, band_group, "PIXEL_SIZE", band_string, 1,
        &band_metadata->pixel_size) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the image start time for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->image_start_time),
        IAS_ODL_Double, odl_data, band_group, "START_TIME", band_string, 1,
        &band_metadata->image_start_time) != SUCCESS)
    {
        return ERROR;
    }

    /* Read the line time for this band */
    if (construct_and_get_odl_field(sizeof(band_metadata->seconds_per_line),
        IAS_ODL_Double, odl_data, band_group, "LINE_TIME", band_string, 1,
        &band_metadata->seconds_per_line) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite mean height */
    if (construct_and_get_odl_field(sizeof(
        band_metadata->satellite.mean_height), IAS_ODL_Double, odl_data,
        band_group, "MEAN_HEIGHT", band_string, 1, 
        &band_metadata->satellite.mean_height) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite mean L1r line/samp offsets */
    if (construct_and_get_odl_field(sizeof(offsets), IAS_ODL_Double, odl_data, 
        band_group, "MEAN_L1R_LINE_SAMP", band_string, 2, offsets) != SUCCESS)
    {
        return ERROR;
    }
    band_metadata->satellite.line_terms.l1r_mean_offset = offsets[0];
    band_metadata->satellite.samp_terms.l1r_mean_offset = offsets[1];

    /* Satellite mean L1T line/samp offsets */
    if (construct_and_get_odl_field(sizeof(offsets), IAS_ODL_Double, odl_data,
        band_group, "MEAN_L1T_LINE_SAMP", band_string, 2, offsets) != SUCCESS)
    {
        return ERROR;
    }
    band_metadata->satellite.line_terms.l1t_mean_offset = offsets[0];
    band_metadata->satellite.samp_terms.l1t_mean_offset = offsets[1];

    /* Satellite mean vector */
    if (construct_and_get_odl_field(sizeof(
        band_metadata->satellite.mean_offset), IAS_ODL_Double, odl_data, 
        band_group, "MEAN_SAT_VECTOR", band_string, 3, 
        &band_metadata->satellite.mean_offset) != SUCCESS)
    {
        return ERROR;
    }

    /* L1T image corner lines */
    if (construct_and_get_odl_field(sizeof(
        band_metadata->active_l1t_corner_lines), IAS_ODL_Double, 
        odl_data, band_group, "L1T_IMAGE_CORNER_LINES", band_string, 4,
        &band_metadata->active_l1t_corner_lines) != SUCCESS)
    {
        return ERROR;
    }

    /* L1T image corner samples */
    if (construct_and_get_odl_field(sizeof(
        band_metadata->active_l1t_corner_samps), IAS_ODL_Double, 
        odl_data, band_group, "L1T_IMAGE_CORNER_SAMPS", band_string, 4,
        &band_metadata->active_l1t_corner_samps) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite x numerator coefficients */
    if (construct_and_get_odl_field(num_size, IAS_ODL_Double, odl_data,
        band_group, "SAT_X_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF, 
        band_metadata->satellite.x_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite x denominator coefficients */
    if (construct_and_get_odl_field(den_size, IAS_ODL_Double, odl_data,
        band_group, "SAT_X_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->satellite.x_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite y numerator coefficients */
    if (construct_and_get_odl_field(num_size, IAS_ODL_Double, odl_data,
        band_group, "SAT_Y_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->satellite.y_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite y denominator coefficients */
    if (construct_and_get_odl_field(den_size, IAS_ODL_Double, odl_data,
        band_group, "SAT_Y_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->satellite.y_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite z numerator coefficients */
    if (construct_and_get_odl_field(num_size, IAS_ODL_Double, odl_data,
        band_group, "SAT_Z_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->satellite.z_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Satellite z denominator coefficients */
    if (construct_and_get_odl_field(den_size, IAS_ODL_Double, odl_data,
        band_group, "SAT_Z_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->satellite.z_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar vector RPC model */
    band_metadata->solar.mean_height = band_metadata->satellite.mean_height;
    band_metadata->solar.line_terms.l1r_mean_offset 
        = band_metadata->satellite.line_terms.l1r_mean_offset;
    band_metadata->solar.samp_terms.l1r_mean_offset 
        = band_metadata->satellite.samp_terms.l1r_mean_offset;
    band_metadata->solar.line_terms.l1t_mean_offset 
        = band_metadata->satellite.line_terms.l1t_mean_offset;
    band_metadata->solar.samp_terms.l1t_mean_offset 
        = band_metadata->satellite.samp_terms.l1t_mean_offset;

    /* Solar mean vector */
    if (construct_and_get_odl_field(sizeof(band_metadata->solar.mean_offset),
        IAS_ODL_Double, odl_data, band_group, "MEAN_SUN_VECTOR", band_string, 3,
        &band_metadata->solar.mean_offset) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar x numerator coefficients */
    if (construct_and_get_odl_field(num_size, IAS_ODL_Double, odl_data,
        band_group, "SUN_X_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->solar.x_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar x denominator coefficients */
    if (construct_and_get_odl_field(den_size, IAS_ODL_Double, odl_data,
        band_group, "SUN_X_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->solar.x_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar y numerator coefficients */
    if (construct_and_get_odl_field(num_size, IAS_ODL_Double, odl_data,
        band_group, "SUN_Y_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->solar.y_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar y denominator coefficients */
    if (construct_and_get_odl_field(den_size, IAS_ODL_Double, odl_data,
        band_group, "SUN_Y_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->solar.y_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar z numerator coefficients */
    if (construct_and_get_odl_field(num_size, IAS_ODL_Double, odl_data,
        band_group, "SUN_Z_NUM_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF,
        band_metadata->solar.z_terms.numerator) != SUCCESS)
    {
        return ERROR;
    }

    /* Solar z denominator coefficients */
    if (construct_and_get_odl_field(den_size, IAS_ODL_Double, odl_data,
        band_group, "SUN_Z_DEN_COEF", band_string, IAS_ANGLE_GEN_ANG_RPC_COEF 
        - 1, band_metadata->solar.z_terms.denominator) != SUCCESS)
    {
        return ERROR;
    }

    /* SCA List */
    if (construct_and_get_odl_field(IAS_MAX_NSCAS * sizeof(int), 
        IAS_ODL_Int, odl_data, band_group, "SCA_LIST", band_string, 
        band_metadata->num_scas, sca_list) != SUCCESS)
    {
        return ERROR;
    }

    /* Loop through the SCAs */
    for (sca = 0; sca < band_metadata->num_scas; sca++)
    {
        if (ias_angle_gen_read_ang_sca(odl_data, band_metadata->band_number,
            sca_list[sca], band_group, &band_metadata->sca_metadata[sca]) 
            != SUCCESS)
        {
            IAS_LOG_ERROR("Reading the SCA metadata for SCA %d", sca);
            return ERROR;
        }
    }

    return SUCCESS;
}

/*******************************************************************************
Name: ias_angle_gen_read_ang

Purpose: Read the ANG metadata file.

Returns: 
    Type = integer
    SUCCESS / ERROR
 ******************************************************************************/
int ias_angle_gen_read_ang
(
    const char *ang_filename,        /* I: Angle file name to read */
    IAS_ANGLE_GEN_METADATA *metadata /* O: Metadata structure to load */
)       
{
    IAS_OBJ_DESC *odl_data; /* Metadata ODL object */
    int index;              /* Loop index */

    /* Open file */
    odl_data = ias_odl_read_tree(ang_filename);
    if (!odl_data)
    {
        IAS_LOG_ERROR("Opening input metadata file %s", ang_filename);
        return ERROR;
    }

    /* Read data */
    /* File information group */
    if (ias_angle_gen_read_ang_header(odl_data, metadata) != SUCCESS)
    {
        IAS_LOG_ERROR("Reading the file header from %s", ang_filename);
        ias_odl_free_tree(odl_data);
        return ERROR;
    }

    /* Projection information group */
    if (ias_angle_gen_read_ang_projection(odl_data, metadata) 
        != SUCCESS)
    {
        IAS_LOG_ERROR("Reading the projection information from %s", 
            ang_filename);
        ias_odl_free_tree(odl_data);
        return ERROR;
    }

    /* Ephemeris data group and solar vector data group */
    if (ias_angle_gen_read_ang_ephemeris(odl_data, metadata) != SUCCESS)
    {
        IAS_LOG_ERROR("Reading ephemeris data from %s", ang_filename);
        ias_odl_free_tree(odl_data);
        return ERROR;
    }

    /* Band information */
    for (index = 0; index < IAS_MAX_NBANDS; index++)
    {
        if (!metadata->band_present[index]) continue;

        if (ias_angle_gen_read_ang_band(odl_data, 
            &metadata->band_metadata[index]) != SUCCESS)
        {
            IAS_LOG_ERROR("Reading band metadata for index %d from %s", index, 
                ang_filename);
            ias_odl_free_tree(odl_data);
            ias_angle_gen_free(metadata);
            return ERROR;
        }
    }

    /* Release the ODL structure */
    ias_odl_free_tree(odl_data);

    return SUCCESS;
}
