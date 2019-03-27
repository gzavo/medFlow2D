#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "utils/interpolation.h"
#include "utils/miniz.h"
#include "ini/INIReader.h"
//#include "utils/tinycthread.h"

#define NUMPREC 15		// digit precision for text file save/load

#define EPS_SLOW 1e-9

typedef double real;

using namespace std;

int numOfCPUs = 2;

string inBaseDir;
string inBaseFileName;
string outBaseDir;
string outBaseFileName;
string outletFileName;

string startPositionsFile;

int numOfPoints = 0;

int idxFrom = 0;
int idxTo = 0;
int idxStep = 0;
int numTimeSlices = 0;

int geometryX = 0;
int geometryY = 0;

real cfdDt = 0.0;
real cfdDx = 0.0;

real tracerDt = 0;
real dtBetweenTimeSlices = 0.0;
real fullTracingTime = 0;

int currentHead = 0;
real currentTimeBetweenSnapshots = 0.0;
real currentTimeRatioBetweenSnapshots = 0;
real currentTime = 0;

real saveEvery = 1000;

typedef struct {
    real x,y;
} CPoint;

typedef struct {
    CPoint p1, p2;
} Outlet;

typedef struct {
    int fromIdx;
    int toIdx;
} WorkerData;

int numOfOutlets = 0;
Outlet *out; // Outlet boundaries

CPoint *p;  // Tracer particle coordinates
CPoint **v; // Velocity snapshots from medFlow2D
real *residence_time;
int *finalOutlet;

inline bool fileExists (const std::string& name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

string zeroPadNumber(int num)
{
    std::ostringstream ss;
    ss << std::setw( 7 ) << std::setfill( '0' ) << num;
    return ss.str();
}

bool readSetupIni(const string iniFileName)
{
    cout << "Loading ini file: " << iniFileName << endl;

    if (!fileExists(iniFileName)) {
        std::cout << "Cannot open file: " << iniFileName << "!\n";
        return false;
    }

    INIReader ini(iniFileName);

    if (ini.ParseError() < 0) {
        std::cout << "Cannot parse " << iniFileName << " as a valid INI file!\n";
        return false;
    }

    // Read values from the ini file

    // cfd
    inBaseDir = ini.Get("cfd", "inBaseDir", "./");
    inBaseFileName = ini.Get("cfd", "inBaseFileName", "full_");
    idxFrom = ini.GetInteger("cfd", "idxFrom", 0);
    idxTo = ini.GetInteger("cfd", "idxTo", 0);
    idxStep = ini.GetInteger("cfd", "idxStep", 0);
    cfdDt = ini.GetReal("cfd", "cfdDt", 0.0);
    cfdDx = ini.GetReal("cfd", "cfdDx", 0.0);


    // tracer
    startPositionsFile = ini.Get("tracer", "pointListFile", "./points.txt");
    numOfPoints = ini.GetInteger("tracer", "numOfPoints", 0);
    tracerDt = ini.GetReal("tracer", "tracerDt", 0.0001);
    fullTracingTime = ini.GetReal("tracer", "fullTracingTime", 0.0);
    outletFileName = ini.Get("tracer", "outletsFile", "./outlets.txt");
    numOfOutlets = ini.GetInteger("tracer", "numOfOutlets", 0);

    // output
    outBaseDir = ini.Get("output", "outBaseDir", "./");
    outBaseFileName = ini.Get("output", "outBaseFileName", "res_");
    saveEvery = ini.GetReal("output", "saveEvery", 10000.0);

    return true;
}

bool readStartingPositions(const string fileName)
{
    cout << "Loading positions file: " << fileName << endl;

    if (!fileExists(fileName)) {
        std::cout << "Cannot open file: " << fileName << "!\n";
        return false;
    }

    p = new CPoint[numOfPoints];

    ifstream fpos(fileName);

    for(int i = 0; i < numOfPoints; i++)
        fpos >> p[i].x >> p[i].y;

    return true;
}

bool readOutlets(const string fileName)
{
    cout << "Loading outlets file: " << fileName << endl;

    if (!fileExists(fileName)) {
        std::cout << "Cannot open file: " << fileName << "!\n";
        return false;
    }

    out = new Outlet[numOfOutlets];

    ifstream fpos(fileName);

    for(int i = 0; i < numOfOutlets; i++)
        fpos >> out[i].p1.x >> out[i].p1.y >> out[i].p2.x >>out[i].p2.y;

    return true;

}

CPoint *readCompressedVelocityfield(const string fileName)
{
    if (!fileExists(fileName)) {
        std::cout << "Cannot open file: " << fileName << "!\n";
        return NULL;
    }

    mz_bool status;
    mz_zip_archive zip_archive;

    size_t uncomp_size;
    char *cData;

    // Extracting txt file name

    string txtFileName = fileName;

    // Remove directory if present.
    // Do this before extension removal in case the directory has a period character.
    const size_t last_slash_idx = txtFileName.find_last_of("\\/");
    if (std::string::npos != last_slash_idx)
    {
        txtFileName.erase(0, last_slash_idx + 1);
    }

    // Remove extension if present. (.zip)
    const size_t period_idx = txtFileName.rfind('.');
    if (std::string::npos != period_idx)
    {
        txtFileName.erase(period_idx);
    }

    cout << "Text archive name: " << txtFileName << endl;

    // Reading compressed file
    memset(&zip_archive, 0, sizeof(zip_archive));

    status = mz_zip_reader_init_file(&zip_archive, fileName.c_str(), 0);
    if (!status)
    {
        cout << "Reading compressed file " << fileName << " failed!" << endl;
        return NULL;
    }

    cData = (char *)mz_zip_reader_extract_file_to_heap(&zip_archive, txtFileName.c_str(), &uncomp_size, 0);
    if (!cData)
    {
        printf("mz_zip_reader_extract_file_to_heap() failed!\n");
        mz_zip_reader_end(&zip_archive);
        return NULL;
    }

    // Process currentVelocityField (also, change it to pixel/second)
    string input(cData);
    istringstream inStream( input );

    int resX, resY;
    inStream >> resX >> resY;

    if(geometryX == 0 && geometryY == 0)
    {
        geometryX = resX;
        geometryY = resY;
    }
    else
        if(geometryX != resX || geometryY != resY)
        {
            cout << "The resolutions of the currentVelocityField files do not match!" << endl;
            return NULL;
        }

    cout << "Geometry: " << resX << ", " << resY << endl;

    CPoint *currentVelocityField = new CPoint[resX*resY];

    for(int i = 0; i < resX*resY; i++) {
        real t1, t2, t3, t4, t5;
        inStream >> currentVelocityField[i].x >> currentVelocityField[i].y >> t1 >> t2 >> t3 >> t4 >> t5;
        currentVelocityField[i].x /= cfdDx; currentVelocityField[i].y /= cfdDx;     // scale to pixel/s units
    }

    mz_zip_reader_end(&zip_archive);


    return currentVelocityField;
}

inline double velMagSqr(const CPoint &p)
{
    return p.x*p.x + p.y*p.y;
}

inline double velMag(const CPoint &p)
{
    return sqrt(velMagSqr(p));
}

inline int gT(const int x, const int y)
{
    return geometryX*y+x;
}


inline void interpSpaceCub(CPoint *resVel, const CPoint *velocityField, const CPoint *x) {
    int xm, ym, zm, xp, yp, zp;

    // Get surrounding corner index positions
    xm = (int) floor(x->x);
    ym = (int) floor(x->y);


    real dist_x = (x->x - xm);
    real dist_y = (x->y - ym);

    real w[16];

    interpolation_cubic_2D_weights(dist_x, dist_y, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7], w[8], w[9], w[10], w[11], w[12], w[13], w[14], w[15] );
    resVel->x = 0; resVel->y = 0;

    int scount = 0;
    for(int iy= ym-1; iy< ym+3; iy++)
        for(int ix= xm-1; ix< xm+3; ix++) {
            int idx = gT(ix, iy);
            resVel->x += velocityField[idx].x*w[scount];
            resVel->y += velocityField[idx].y*w[scount];
            scount++;
        }
}

inline void interpSpaceLin(CPoint *resVel, const CPoint *velocityField, const CPoint *x) {

    real xm, ym, xp, yp;

    // Get surrounding corner index positions
    xm = floor(x->x);
    xp = xm + 1;
    ym = floor(x->y);
    yp = ym + 1;

    CPoint corn[] = {    {xm,ym},
                         {xp,ym},
                         {xm,yp},
                         {xp,yp} };
    real w[4];

    real dist_x = (x->x - xm);
    real dist_y = (x->y - ym);

    interpolation_lerp_2D_weights(dist_x, dist_y, w[0], w[1], w[2], w[3]);

    resVel->x = 0; resVel->y = 0;


    // Normalized vectorial sum of corner velocities
    for(int i=0; i < 4; i++) {
        int idx = gT((int)corn[i].x, (int)corn[i].y);
        resVel->x += velocityField[idx].x*w[i];
        resVel->y += velocityField[idx].y*w[i];
    }

}

void stepTimeForward(real dt)
{
    if(currentTimeBetweenSnapshots + dt > dtBetweenTimeSlices)
    {
        currentTimeBetweenSnapshots = fmod(currentTimeBetweenSnapshots + dt, dtBetweenTimeSlices);
        currentHead += 1;
        if (currentHead == numTimeSlices-1)
            currentHead =0;

    }
    else
        currentTimeBetweenSnapshots += dt;

    currentTimeRatioBetweenSnapshots = currentTimeBetweenSnapshots / dtBetweenTimeSlices;
    currentTime += dt;
}

void stepTimeBackward(real dt)
{
    if(currentTimeBetweenSnapshots - dt < 0.0)
    {
        currentTimeBetweenSnapshots = dtBetweenTimeSlices + (currentTimeBetweenSnapshots - dt);
        currentHead -= 1;
        if (currentHead < 0)
            currentHead = numTimeSlices - 2;

    }
    else
        currentTimeBetweenSnapshots -= dt;

    currentTimeRatioBetweenSnapshots = currentTimeBetweenSnapshots / dtBetweenTimeSlices;
    currentTime -= dt;
}

void stepTimeForwardTemp(real dt, int *currentHead, real *currentTimeBetweenSnapshots, real * currentTimeRatioBetweenSnapshots)
{
    if(*currentTimeBetweenSnapshots + dt > dtBetweenTimeSlices)
    {
        *currentTimeBetweenSnapshots = fmod(*currentTimeBetweenSnapshots + dt, dtBetweenTimeSlices);
        currentHead += 1;
        if (*currentHead == numTimeSlices-1)
            *currentHead = 0;

    }
    else
        *currentTimeBetweenSnapshots += dt;

    *currentTimeRatioBetweenSnapshots = *currentTimeBetweenSnapshots / dtBetweenTimeSlices;
    //currentTime += dt;
}

void stepTimeBackwardTemp(real dt, int *currentHead, real *currentTimeBetweenSnapshots, real * currentTimeRatioBetweenSnapshots)
{
    if(*currentTimeBetweenSnapshots - dt < 0.0)
    {
        *currentTimeBetweenSnapshots = dtBetweenTimeSlices + (*currentTimeBetweenSnapshots - dt);
        currentHead -= 1;
        if (*currentHead < 0)
            *currentHead = numTimeSlices - 2;

    }
    else
        *currentTimeBetweenSnapshots -= dt;

    *currentTimeRatioBetweenSnapshots = *currentTimeBetweenSnapshots / dtBetweenTimeSlices;
    //currentTime -= dt;
}


void inline evalVelocity(const CPoint &x, CPoint *resVel, const int currentHead, const real currentTimeRatioBetweenSnapshots) {

    CPoint v1, v2;

    interpSpaceLin(&v1, v[currentHead], &x);
    interpSpaceLin(&v2, v[currentHead+1], &x);

    resVel->x = interpolation_lerp_1D(v1.x, v2.x, currentTimeRatioBetweenSnapshots);
    resVel->y = interpolation_lerp_1D(v1.y, v2.y, currentTimeRatioBetweenSnapshots);
}

inline int passedOutlet(const CPoint &x)
{
    for(int i = 0; i < numOfOutlets; i++)
    {
        if(x.x >= out[i].p1.x && x.x <= out[i].p2.x)
            if(x.y >= out[i].p1.y && x.y <= out[i].p2.y)
                return i;
    }

    return -1;
}

// Returns true if the integrated point is SLOW
bool integrateEuler(const real dt, CPoint &p)
{
    CPoint vel = { 0, 0 };

    evalVelocity(p, &vel, currentHead, currentTimeRatioBetweenSnapshots);

    p.x += vel.x*dt;
    p.y += vel.y*dt;

    if (sqrt(velMagSqr(vel)) < EPS_SLOW)
        return true;

    return false;
}

bool integrateRK4(const real dt, CPoint &p)
{
    //TODO: check for out of bound error at every step > 1

    //Creating local copy of the time state description
    int t_currentHead = currentHead;
    real t_currentTimeBetweenSnapshots = currentTimeBetweenSnapshots;
    real t_currentTimeRatioBetweenSnapshots = currentTimeRatioBetweenSnapshots;


    // first step
    CPoint dp = { 0, 0 };
    evalVelocity(p, &dp, t_currentHead, t_currentTimeRatioBetweenSnapshots);
    //dp is just v yet
    bool isSlow = (sqrt(velMagSqr(dp)) < EPS_SLOW);
    dp.x *= tracerDt; dp.y *= tracerDt;

    //second step
    CPoint dp2 = { 0, 0 }; CPoint p2 = { p.x + dp.x/2.0, p.y + dp.y/2.0 };
    stepTimeForwardTemp(tracerDt / 2.0, &t_currentHead, &t_currentTimeBetweenSnapshots, &t_currentTimeRatioBetweenSnapshots);
    evalVelocity(p2, &dp2,t_currentHead, t_currentTimeRatioBetweenSnapshots);
    dp2.x *= tracerDt; dp2.y *= tracerDt;

    //third step
    CPoint dp3 = { 0, 0 }; CPoint p3 = { p.x + dp2.x / 2.0, p.y + dp2.y / 2.0};
    evalVelocity(p3, &dp3, t_currentHead, t_currentTimeRatioBetweenSnapshots);
    dp3.x *= tracerDt; dp3.y *= tracerDt;

    //fourth step
    CPoint dp4 = { 0, 0 }; CPoint p4 = { p.x + dp3.x, p.y + dp3.y };
    stepTimeForwardTemp(tracerDt / 2.0, &t_currentHead, &t_currentTimeBetweenSnapshots, &t_currentTimeRatioBetweenSnapshots);
    evalVelocity(p4, &dp4, t_currentHead, t_currentTimeRatioBetweenSnapshots);
    dp4.x *= tracerDt; dp4.y *= tracerDt;

    // Move the clock back
    //stepTimeBackward(tracerDt);

    p.x += (dp.x + 2 * dp2.x + 2 * dp3.x + dp4.x) / 6.0;
    p.y += (dp.y + 2 * dp2.y + 2 * dp3.y + dp4.y) / 6.0;

    return isSlow;
}

bool savePositions(const string fileName)
{
    ofstream ofile(fileName);

    for(int i=0; i < numOfPoints; i++)
    {
        ofile << p[i].x << " "  << p[i].y << " " << endl;
    }

    ofile.close();

    return true;
}

bool saveResidence(const string fileName)
{
    ofstream ofile(fileName);

    for(int i=0; i < numOfPoints; i++)
    {
        ofile << residence_time[i] << " "  << finalOutlet[i] << endl;
    }

    ofile.close();

    return true;
}

int integratorWorker(void *arg) {

    WorkerData *w = (WorkerData*)arg;

    for(int i = w->fromIdx; i < w->toIdx; i++) {
        if (residence_time[i] == 0) {
            bool isSlow = integrateRK4(tracerDt, p[i]);

            if (passedOutlet(p[i]) > -1) {
                finalOutlet[i] = passedOutlet(p[i]);
                residence_time[i] = currentTime;
            }
            else if (isSlow)
                residence_time[i] = -1.0;
        }
    }
}

int main(int argc, char *argv[]) {
    // Verify input arguments
    if (argc != 2) {
        std::cout << "Usage: " << argv[0]
        << " setup.ini" << std::endl;
        return EXIT_FAILURE;
    }

    // Read ini file
    if(!readSetupIni(argv[1]))
        return -1;

    // Calculating necessary values
    numTimeSlices = round((idxTo - idxFrom)/(real)idxStep)+1;
    cout << "Number of time slices: "  << numTimeSlices << endl;

    dtBetweenTimeSlices = idxStep*cfdDt;
    cout << "Input flow field of " << dtBetweenTimeSlices * (numTimeSlices-1) << "s." << endl;

    // Read starting positions
    if(!readStartingPositions(startPositionsFile))
        return -1;

    // Read outlets
    if(!readOutlets(outletFileName))
        return -1;

    // Read compressed velocity field
    v = new CPoint*[numTimeSlices];

    for(int i = 0; i < numTimeSlices; i++)
    {
        string compFileName = inBaseDir + "/" + inBaseFileName + zeroPadNumber(idxFrom + i*idxStep) + ".txt.zip";
        cout << "Loading archive: (" << i+1 << "/" << numTimeSlices << ") " << compFileName << endl;
        v[i] = readCompressedVelocityfield(compFileName);

        if(v[i]==NULL)
        {
            cout << "File loading error, exiting..." << endl;
            return -1;
        }
    }

    residence_time = new real[numOfPoints];
    finalOutlet = new int[numOfPoints];

    for(int i=0; i < numOfPoints; i++)
    {
        residence_time[i] = 0;
        finalOutlet[i] = -1;
    }

    int maxIter = (int)ceil(fullTracingTime / tracerDt);

    int saveIdx = (int)round(saveEvery / tracerDt);

    cout << "Starting trace for " << maxIter << " steps..." << endl;

    currentTime = 0;

//    int workPerProcess = (int)floor((real) numOfPoints / numOfCPUs);
//    WorkerData *pdata = new WorkerData[numOfCPUs];
//    thrd_t *t = new thrd_t[numOfCPUs];

    // main loop
    for(int iT = 0; iT < maxIter+1; iT++)
    {
        // Launching threads
//        for(int i = 0; i < numOfCPUs; i++) {
//
//            pdata[i].fromIdx = i * workPerProcess; //1 is the first cell (see ghost cells)
//
//            if (i == numOfCPUs - 1)  //last thread does some overwork
//                pdata[i].toIdx = numOfPoints;
//            else
//                pdata[i].toIdx = (i + 1) * workPerProcess;
//
//            thrd_create(&t[i], &integratorWorker, (void*)&pdata[i]);
//        }
//
//        // Wait for threads to finish
//        for(int i=0; i < numOfCPUs; i++) {
//            thrd_join(t[i], NULL);
//        }
        #pragma omp parallel for
        for(int i = 0; i < numOfPoints; i++) {
            if (residence_time[i] == 0) {
                bool isSlow = integrateRK4(tracerDt, p[i]);

                if (passedOutlet(p[i]) > -1) {
                    finalOutlet[i] = passedOutlet(p[i]);
                    residence_time[i] = currentTime;
                }
                else if (isSlow)
                    residence_time[i] = -1.0;
            }
        }

        stepTimeForward(tracerDt);

        if(iT%saveIdx == 0)
        {
            std::cout << "Time: " << currentTime << " s / " << fullTracingTime << " s  Iter: " << iT << endl;

            string outFile_pos = outBaseDir + "/" + outBaseFileName + "pos_" + zeroPadNumber(iT) + ".txt";
            string outFile_res = outBaseDir + "/" + outBaseFileName + "res_" + zeroPadNumber(iT) + ".txt";

            savePositions(outFile_pos);
            saveResidence(outFile_res);
        }

    }

    cout << "Done tracing! :)" << endl;

    return 0;
}