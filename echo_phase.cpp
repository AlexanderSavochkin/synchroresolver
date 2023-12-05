/*
Writing wav file is taken from here: https://github.com/Thrifleganger/audio-programming-youtube/blob/master/wav-file-format/main.cpp
*/

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <portaudio.h>
#include "pa_ringbuffer.h"

#define PA_SAMPLE_TYPE      paFloat32


const float SAMPLE_RATE = 44100;
const unsigned int FRAMES_PER_BUFFER = 64;
//This frequency is close to 400HZ, but has integer number of samples per period (110 samples per period)
const float DEFAULT_FREQUENCY_HZ = 400.90909090909093f;
const unsigned int SAMPLES_PER_PERIOD = 110;


const float SMOOTHING_PARAM = 0.01f;
const float PI = 3.1415927f;
const size_t RING_BUFFER_SIZE = 16384;
const size_t RING_BUFFER_CHUNK_SIZE = 1024;
const int CIRCULAR_BUFFER_POLL_PERIOD = 10;

using std::vector;
using std::string;
using std::ofstream;
using std::clog;
using std::endl;
using std::ios;
using std::cos;
using std::sin;
using std::stoi;
using std::min;

void reportDevices()
{
    int numDevices = Pa_GetDeviceCount();
    clog << "Found " << numDevices << " audio devices" << endl;
    if (numDevices > 0)
    {
        for (int i = 0; i < numDevices; ++i)
        {
            const PaDeviceInfo* deviceInfo = Pa_GetDeviceInfo(i);
            clog
                << i << ' '
                << "------" << deviceInfo->name << "-------" << endl
                << "Maximum input channels: " << deviceInfo->maxInputChannels << endl
                << "Maximum output channels: " << deviceInfo->maxOutputChannels << endl
                << "Default low input latency: " << deviceInfo->defaultLowInputLatency << endl
                << "Default high input latency: " << deviceInfo->defaultHighInputLatency << endl
                << endl;
        }
    }
    else
    {
        clog << "PortAudio get device count error: " << Pa_GetErrorText(numDevices) << endl;
	exit(1);
    }
}

struct FeedBackData
{
    float T;

    float CosineSum0;
    float SineSum0;

    float CosineSum1;
    float SineSum1;

    float Frequency;

    int running;
    int numPeriodsToCalculateAngle;

    vector<float> inputRingBufferData;
    PaUtilRingBuffer inputRingBuffer;

    FeedBackData() : 
        T(0.0f),
        CosineSum0(0.0f),
        SineSum0(0.0f),
        CosineSum1(0.0f),
        SineSum1(0.0f),
        Frequency(DEFAULT_FREQUENCY_HZ),
        inputRingBufferData(RING_BUFFER_SIZE),
	numPeriodsToCalculateAngle(20),
        running(0)
	{        
        long initRingBufferResult = PaUtil_InitializeRingBuffer(
            &inputRingBuffer,
            sizeof(float),
            RING_BUFFER_SIZE,
            &inputRingBufferData[0]);

        if (initRingBufferResult != 0) {
            clog << "Buffer size is not power of two:" << RING_BUFFER_SIZE;
            throw std::invalid_argument("Ring buffer size is not power of two");
        }
    }
};

int processAudioDataCallback(
    const void* inputBuffer,
    void* outputBuffer,
    unsigned long framesPerBuffer,
    const PaStreamCallbackTimeInfo* timeInfo,
    PaStreamCallbackFlags statusFlags,
    void* userData)
{
    /* Cast data passed through stream to our structure. */
    FeedBackData* feedBackData = (FeedBackData*)userData;
    float* out = (float*)outputBuffer;
    const float* inp = (const float*)inputBuffer;
    unsigned int i;

    static float dCosineSum0 = 0.0f;
    static float dSineSum0 = 0.0f;
    static float dCosineSum1 = 0.0f;
    static float dSineSum1 = 0.0f;
    static int numFramesAccumulated = 0;

    ring_buffer_size_t elementsWriteable = PaUtil_GetRingBufferWriteAvailable(&feedBackData->inputRingBuffer);

    PaUtil_WriteRingBuffer(
        &feedBackData->inputRingBuffer,
        inputBuffer,
        min((ring_buffer_size_t)(2  /*2 Channels*/ * framesPerBuffer), elementsWriteable));

    for (i = 0; i < framesPerBuffer; i++)
    {
        float outCosine = cos(float(2 * PI * feedBackData->Frequency * feedBackData->T) / SAMPLE_RATE);
        float outSine = sin(float(2 * PI * feedBackData->Frequency * feedBackData->T) / SAMPLE_RATE);

        *out++ = outCosine;  /* left */
        *out++ = outSine;   /* right */

        dCosineSum0 += (outCosine * (*inp));
        dSineSum0 += (outSine * (*inp));
        ++inp;

        dCosineSum1 += (outCosine * (*inp));
        dSineSum1 += (outSine * (*inp));
        ++inp;

	numFramesAccumulated++;
	if (numFramesAccumulated == feedBackData->numPeriodsToCalculateAngle * SAMPLES_PER_PERIOD)
	{
            feedBackData->CosineSum0 = dCosineSum0;
            feedBackData->SineSum0 = dSineSum0;
            feedBackData->CosineSum1 = dCosineSum1;
            feedBackData->SineSum1 = dSineSum1;

	    numFramesAccumulated = 0;

	    dCosineSum0 = 0.0f;
	    dSineSum0 = 0.0f;
	    dCosineSum1 = 0.0f;
	    dSineSum1 = 0.0f;
	}

	feedBackData->T += 1.0f / SAMPLE_RATE;

        if (feedBackData->T > (SAMPLE_RATE / feedBackData->Frequency))
        {
            feedBackData->T -= (SAMPLE_RATE / feedBackData->Frequency);
        }
    }


    feedBackData->running++;

    return paContinue;
}

FeedBackData feedBackData;

void initializePortAudio()
{
    PaError err = Pa_Initialize();
    if (err != paNoError)
    {
        clog << "PortAudio initialization error: " << Pa_GetErrorText(err) << endl;
        exit(1);
    }
}

void terminatePortAudio()
{
    PaError err = Pa_Terminate();
    if (err != paNoError)
    {
        clog << "PortAudio termination error: " << Pa_GetErrorText(err) << endl;
        exit(1);
    }
}

void openStream(PaStream** stream, const PaStreamParameters* inputParameters, const PaStreamParameters* outputParameters, FeedBackData* feedBackData)
{
    PaError err = Pa_OpenStream(
        stream,
        inputParameters,
        outputParameters,
        SAMPLE_RATE,
        FRAMES_PER_BUFFER,
        paNoFlag, //flags that can be used to define dither, clip settings and more
        processAudioDataCallback, //your callback function
        feedBackData);

    if (err != paNoError)
    {
        clog << "Error opening PortAudio stream: " << Pa_GetErrorText(err) << endl;
        exit(1);
    }
}

void closeStream(PaStream* stream)
{
    PaError err = Pa_CloseStream(stream);
    if (err != paNoError)
    {
        clog << "Error closing PortAudio stream: " << Pa_GetErrorText(err) << endl;
        exit(1);
    }
}

void startStream(PaStream* stream)
{
    PaError err = Pa_StartStream(stream);
    if (err != paNoError)
    {
        clog << "Error starting PortAudio stream: " << Pa_GetErrorText(err) << endl;
        exit(1);
    }
}

void stopStream(PaStream* stream)
{
    clog << "Stopping PortAudio stream" << endl;
    PaError err = Pa_StopStream(stream);
    if (err != paNoError)
    {
        clog << "Error stopping PortAudio stream: " << Pa_GetErrorText(err) << endl;
        exit(1);
    }
    clog << "Stopped PortAudio stream" << endl;
}

void sleep(int milliseconds)
{
    Pa_Sleep(milliseconds);
}

void printFeedbackData(const FeedBackData& feedBackData)
{
    auto dotprod = feedBackData.CosineSum0 * feedBackData.CosineSum1 + feedBackData.SineSum0 * feedBackData.SineSum1;
    auto amplitude0 = sqrt(feedBackData.CosineSum0 * feedBackData.CosineSum0 + feedBackData.SineSum0 * feedBackData.SineSum0);
    auto amplitude1 = sqrt(feedBackData.CosineSum1 * feedBackData.CosineSum1 + feedBackData.SineSum1 * feedBackData.SineSum1);
    //auto angle = dotprod / (amplitude0 * amplitude1);
    auto angle = atan2(feedBackData.SineSum0, feedBackData.CosineSum0) - atan2(feedBackData.SineSum1, feedBackData.CosineSum1);
    clog << "alpha = " << angle << " Cos0 = " << feedBackData.CosineSum0 << " Sin0 = " << feedBackData.SineSum0 << " Cos1 = " << feedBackData.CosineSum1 << " Sin1 = " << feedBackData.SineSum1 << " Amplitute0 = " << amplitude0 << " Amplitude1 = " << amplitude1 << ' ' << feedBackData.running << "                                          ";
}

void printUsage()
{
    clog << "Usage: program [options]" << endl;
    clog << "Options:" << endl;
    clog << "  -l              List available audio devices" << endl;
    clog << "  -i <deviceID>   Specify the input device ID" << endl;
    clog << "  -o <deviceID>   Specify the output device ID" << endl;
    clog << "  -f <frequency>  Set the frequency in Hz (default: 400)" << endl;
    clog << "  -w <filename>   Write the output to a wav file <filename>" << endl;
    clog << "  -t <time Ms>    Recording time duration <Ms>" << endl;
    clog << "  -h              Display this help message" << endl;
}

void writeToFile(
    ofstream &file,
    int value,
    int size) 
{
    file.write(
        reinterpret_cast<const char*> (&value),
        size);
}

int main(int argc, char* argv[])
{
    PaStream* stream;

    int inputDeviceID = -1;
    int outputDeviceID = -1;
    int frequency = DEFAULT_FREQUENCY_HZ;

    string wavFileName = "";
    long timeToRecordMs = -1;


    if (argc == 1) 
    {
        printUsage();
        return 0;
    }

    // Process command-line arguments
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-l")
        {
            initializePortAudio();
            reportDevices();
            return 0;
        }
        else if (arg == "-i" && i + 1 < argc)
        {
            inputDeviceID = stoi(argv[i + 1]);
            ++i;
        }
        else if (arg == "-o" && i + 1 < argc)
        {
            outputDeviceID = stoi(argv[i + 1]);
            ++i;
        }
        else if (arg == "-f" && i + 1 < argc)
        {
            frequency = stoi(argv[i + 1]);
            ++i;
        }
        else if (arg == "-w" && i + 1 < argc)
        {
            wavFileName = argv[i + 1];
            ++i;
        }
        else if (arg == "-t" && i + 1 < argc)
        {
            timeToRecordMs = stoi(argv[i + 1]);
	    ++i;
        }
        else if (arg == "-h")
        {
            printUsage();
            return 0;
        }
        else
        {
            clog << "Invalid argument: " << arg << endl;
            printUsage();
            return 1;
        }
    }

    clog << "Input device ID: " << inputDeviceID << endl;
    clog << "Output device ID: " << outputDeviceID << endl;

    initializePortAudio();

    if (inputDeviceID == -1)
    {
        inputDeviceID = Pa_GetDefaultInputDevice();
    }

    if (inputDeviceID == paNoDevice)
    {
        clog << "Error: No default input device." << endl;
        terminatePortAudio();
        return 1;
    }

    PaStreamParameters inputParameters;
    inputParameters.device = inputDeviceID;
    inputParameters.channelCount = 2;       /* stereo input */
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = nullptr;

    if (outputDeviceID == -1)
    {
        outputDeviceID = Pa_GetDefaultOutputDevice();
    }

    if (outputDeviceID == paNoDevice)
    {
        clog << "Error: No default output device." << endl;
        terminatePortAudio();
        return 1;
    }

    PaStreamParameters outputParameters;
    outputParameters.device = outputDeviceID;
    outputParameters.channelCount = 2;       /* stereo output */
    outputParameters.sampleFormat = PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = nullptr;

    feedBackData.Frequency = frequency;

    openStream(&stream, &inputParameters, &outputParameters, &feedBackData);
    startStream(stream);

    {
        ofstream audioFile;
        int preAudioPosition = 0;
        if (!wavFileName.empty())
        {
            audioFile.open(wavFileName, ios::binary);
            if (!audioFile.is_open())
            {
                clog << "Error opening file " << wavFileName << endl;
                return 1;
            }

            audioFile << "RIFF";
            audioFile << "----";
            audioFile << "WAVE";

            // Format chunk
            audioFile << "fmt ";
            writeToFile(audioFile, 16, 4); // Size
            writeToFile(audioFile, 1, 2); // Compression code [format]
            writeToFile(audioFile, 2, 2); // Number of channels
            writeToFile(audioFile, SAMPLE_RATE, 4); // Sample rate
            writeToFile(audioFile, SAMPLE_RATE * 2, 4 ); // Byte rate: 2 byte per sample
            writeToFile(audioFile, /*bitDepth / 8*/ 2, 2); // Block align
            writeToFile(audioFile, /*bitDepth*/ 16, 2); // Bit depth

            //Data chunk
            audioFile << "data";
            audioFile << "----";

            preAudioPosition = audioFile.tellp();
        }

        float chunk_buffer[2 * RING_BUFFER_CHUNK_SIZE];
        unsigned long timestep = 0;	
        while (Pa_IsStreamActive(stream))
        {
            int dataready;
            while((dataready = PaUtil_GetRingBufferReadAvailable(&feedBackData.inputRingBuffer)) < 2 * RING_BUFFER_CHUNK_SIZE)
            {
                //clog << "Data ready: " << dataready << endl;

                sleep(CIRCULAR_BUFFER_POLL_PERIOD);
            }
	    //clog << "Got enough data! " << dataready << endl;

            PaUtil_ReadRingBuffer(
                &feedBackData.inputRingBuffer,
                chunk_buffer,
                2 * RING_BUFFER_CHUNK_SIZE
            );

            if (audioFile.is_open())
            {
                for(int i = 0; i < RING_BUFFER_CHUNK_SIZE; i++ )
                {
                    auto sample = chunk_buffer[2 * i]; 
                    int intSample = static_cast<int> (32000 * sample );
                    writeToFile(audioFile, intSample, 2);

	            sample = chunk_buffer[2 * i + 1]; 
                    intSample = static_cast<int> (32000 * sample);
                    writeToFile(audioFile, intSample, 2);
                }

                int postAudioPosition = audioFile.tellp();

                //Update sizes in the header
                audioFile.seekp(preAudioPosition - 4);
                writeToFile(audioFile, postAudioPosition - preAudioPosition, 4);

                audioFile.seekp(4, ios::beg);
                writeToFile(audioFile, postAudioPosition - 8, 4);

		audioFile.seekp(0, ios::end);
            }

            ++timestep;

            if (timestep % 10 == 0)
            {
                clog << "\r";
                printFeedbackData(feedBackData);
            }

	    if ( timeToRecordMs != -1 && timeToRecordMs < timestep * CIRCULAR_BUFFER_POLL_PERIOD )
	    {
                break;
	    }
        }

	if (audioFile.is_open())
	{
	    audioFile.close();
	}
    }

    stopStream(stream);
    closeStream(stream);
    terminatePortAudio();
    return 0;
}
