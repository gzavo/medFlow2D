/*
 * The MIT License (MIT)
 *  
 *  Copyleft (C) 2015  Sun Dro (a.k.a. kala13x)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in all
 *  copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <time.h>
#include "slog.h"

/* Supported colors */
#define CLR_NORM     "\x1B[0m"
#define CLR_RED      "\x1B[31m"
#define CLR_GREEN    "\x1B[32m"
#define CLR_YELLOW   "\x1B[33m"
#define CLR_BLUE     "\x1B[34m"
#define CLR_NAGENTA  "\x1B[35m"
#define CLR_CYAN     "\x1B[36m"
#define CLR_WHITE    "\x1B[37m"
#define CLR_RESET    "\033[0m"

/* Max size of string */
#define MAXMSG 8196

/* Flags */
static slog_flags slg;


/*
 * get_system_date - Intialize date with system date.
 * Argument is pointer of SystemDate structure.
 */
void get_system_date(SystemDate *mdate)
{
    time_t rawtime;
    struct tm *timeinfo;
    rawtime = time(NULL);
    timeinfo = localtime(&rawtime);

    /* Get System Date */
    mdate->year = timeinfo->tm_year+1900;
    mdate->mon = timeinfo->tm_mon+1;
    mdate->day = timeinfo->tm_mday;
    mdate->hour = timeinfo->tm_hour;
    mdate->min = timeinfo->tm_min;
    mdate->sec = timeinfo->tm_sec;
}


/*
 * Get library version. Function returns version and build number of slog
 * library. Return value is char pointer. Argument min is flag for output
 * format. If min is 0, function returns version in full  format, if flag
 * is 1 function returns only version number, For examle: 1.3.0
 */
const char* slog_version(int min)
{
    static char verstr[128];

    /* Version short */
    if (min) sprintf(verstr, "%d.%d.%d",
        SLOGVERSION_MAX, SLOGVERSION_MIN, SLOGBUILD_NUM);

    /* Version long */
    else sprintf(verstr, "%d.%d build %d (%s)",
        SLOGVERSION_MAX, SLOGVERSION_MIN, SLOGBUILD_NUM, __DATE__);

    return verstr;
}


/*
 * strclr - Colorize string. Function takes color value and string
 * and returns colorized string as char pointer. First argument clr
 * is color value (if it is invalid, function retunrs NULL) and second
 * is string with va_list of arguments which one we want to colorize.
 */
char* strclr(int clr, char* str, ...)
{
    /* String buffers */
    static char output[MAXMSG];
    char string[MAXMSG];

    /* Read args */
    va_list args;
    va_start(args, str);
    vsprintf(string, str, args);
    va_end(args);

    /* Handle colors */
    switch(clr)
    {
        case 0:
            sprintf(output, "%s%s%s", CLR_NORM, string, CLR_RESET);
            break;
        case 1:
            sprintf(output, "%s%s%s", CLR_GREEN, string, CLR_RESET);
            break;
        case 2:
            sprintf(output, "%s%s%s", CLR_RED, string, CLR_RESET);
            break;
        case 3:
            sprintf(output, "%s%s%s", CLR_YELLOW, string, CLR_RESET);
            break;
        case 4:
            sprintf(output, "%s%s%s", CLR_BLUE, string, CLR_RESET);
            break;
        case 5:
            sprintf(output, "%s%s%s", CLR_NAGENTA, string, CLR_RESET);
            break;
        case 6:
            sprintf(output, "%s%s%s", CLR_CYAN, string, CLR_RESET);
            break;
        case 7:
            sprintf(output, "%s%s%s", CLR_WHITE, string, CLR_RESET);
            break;
        default:
            return NULL;
    }

    /* Return output */
    return output;
}


/*
 * log_to_file - Save log in file. Argument aut is string which
 * we want to log. Argument fname is log file path and mdate is
 * SystemDate structure variable, we need it to create filename.
 */
void log_to_file(char *out, const char *fname, SystemDate *mdate)
{
    /* Used variables */
    char filename[PATH_MAX];

    /* Create log filename with date */
    sprintf(filename, "%s-%02d-%02d-%02d.log",
        fname, mdate->year, mdate->mon, mdate->day);

    /* Open file pointer */
    FILE *fp = fopen(filename, "a");
    if (fp == NULL) return;

    /* Write key in file */
    fprintf(fp, "%s", out);

    /* Close file pointer */
    fclose(fp);
}


/*
 * parse_config - Parse config file. Argument cfg_name is path
 * of config file name to be parsed. Function opens config file
 * and parses LOGLEVEL and LOGTOFILE flags from it.
 */
//
//int parse_config(const char *cfg_name)
//{
//    /* Used variables */
//    FILE *file;
//    char *line = NULL;
//    size_t len = 0;
//    ssize_t read;
//    int ret = 1;
//
//    /* Open file pointer */
//    file = fopen(cfg_name, "r");
//    if(file == NULL) return 1;
//
//    /* Line-by-line read cfg file */
//    while ((read = getline(&line, &len, file)) != -1)
//    {
//        /* Find level in file */
//        if(strstr(line, "LOGLEVEL") != NULL)
//        {
//            /* Get log level */
//            slg.level = atoi(line+8);
//            ret = 0;
//        }
//        else if(strstr(line, "LOGTOFILE") != NULL)
//        {
//            /* Get log level */
//            slg.to_file = atoi(line+9);
//            ret = 0;
//        }
//    }
//
//    /* Cleanup */
//    if (line) free(line);
//    fclose(file);
//
//    return ret;
//}


/*
 * Retunr string in slog format. Function takes arguments
 * and returns string in slog format without printing and
 * saveing in file. Return value is char pointer.
 */
char* ret_slog(char *msg, ...)
{
    /* Used variables */
    static char output[MAXMSG];
    char string[MAXMSG];
    SystemDate mdate;

    /* initialise system date */
    get_system_date(&mdate);

    /* Read args */
    va_list args;
    va_start(args, msg);
    vsprintf(string, msg, args);
    va_end(args);

    /* Generate output string with date */
    sprintf(output, "%02d.%02d.%02d-%02d:%02d:%02d - %s",
        mdate.year, mdate.mon, mdate.day, mdate.hour,
        mdate.min, mdate.sec, string);

    /* Return output */
    return output;
}


/*
 * slog - Log exiting process. Function takes arguments and saves
 * log in file if LOGTOFILE flag is enabled from config. Otherwise
 * it just prints log without saveing in file. Argument level is
 * logging level and flag is slog flags defined in slog.h header.
 */
void slog(int level, int flag, const char *msg, ...)
{
    /* Used variables */
    SystemDate mdate;
    char string[MAXMSG];
    char prints[MAXMSG];
    char *output;

    /* Initialise system date */
    get_system_date(&mdate);

    /* Read args */
    va_list args;
    va_start(args, msg);
    vsprintf(string, msg, args);
    va_end(args);

    /* Check logging levels */
    if(level <= slg.level)
    {
        /* Handle flags */
        switch(flag) {
            case 1:
                sprintf(prints, "[LIVE]  %s", string);
                break;
            case 2:
                //sprintf(prints, "[%s]  %s", strclr(1, "INFO"), string);
                sprintf(prints, "[%s]  %s", "INFO", string);
                break;
            case 3:
                //sprintf(prints, "[%s]  %s", strclr(3, "WARN"), string);
                sprintf(prints, "[%s]  %s", "WARN", string);
                break;
            case 4:
                //sprintf(prints, "[%s] %s", strclr(4, "DEBUG"), string);
                sprintf(prints, "[%s] %s", "DEBUG", string);
                break;
            case 5:
                //sprintf(prints, "[%s] %s", strclr(2, "ERROR"), string);
                sprintf(prints, "[%s] %s", "ERROR", string);
                break;
            case 6:
                //sprintf(prints, "[%s] %s", strclr(2, "FATAL"), string);
                sprintf(prints, "[%s] %s", "FATAL", string);
                break;
            case 7:
                //sprintf(prints, "%s", string);
                sprintf(prints, "%s", string);
                break;
            default:
                break;
        }

        /* Print output */
        printf("%s", ret_slog("%s", prints));

        /* Save log in file */
        if (slg.to_file)
        {
            output = ret_slog("%s", string);
            log_to_file(output, slg.fname, &mdate);
        }
    }
}


/*
 * Initialize slog library. Function parses config file and reads log
 * level and save to file flag from config. First argument is file name
 * where log will be saved and second argument conf is config file path
 * to be parsed and third argument lvl is log level for this message.
 */
void init_slog(const char* fname, int useFile, int lvl)
{
    slg.level = lvl;
    slg.fname = fname;
    slg.to_file = useFile;

    /* Parse config file */
    //if (parse_config(conf))
    //{
    //    slog(0, SLOG_WARN, "LOGLEVEL and/or LOGTOFILE flag is not set from config.");

        return;
    //}
}
