#pragma once

#ifndef RECORDING_TOOLS_H
#define RECORDING_TOOLS_H

#include <stdio.h>

#include <queue>
#include <string>

#include <boost/shared_ptr.hpp>

// includes, GL
#include <GL/glew.h>
#include <GL/freeglut.h>
#define BUFFER_OFFSET(bytes) ((GLubyte*) NULL + (bytes))

// define, how much frames are stored in mem before writing to hdd
#define MAX_FRAMES_BEFORE_FLUSH 1

const unsigned char TGAheader[12]={0,0,2,0,0,0,0,0,0,0,0,0};

class ImageBuffer 
{
public:
	union {
		struct { GLint vp[4]; };
		struct { int x, y, width, height; };
	};

	// image data
	GLubyte* data;
	int image_size;

	ImageBuffer(): x(0), y(0), width(0), height(0)
	{};
	ImageBuffer(int _x, int _y, int _width, int _height) : x(_x), y(_y), width(_width), height(_height)
	{};
    ~ImageBuffer()
    {
        delete[] data;
    };

	// initialize
	void init() 
	{ 
		//if( !( vp[0] | vp[1] | vp[2] | vp[3] ) )
		getViewport();
		
		image_size = width * height * 3;
		data = new GLubyte [ image_size ];
	};

	void readPixel() 
	{ 
		glReadPixels(x, y, width, height, GL_BGR, GL_UNSIGNED_BYTE, data);
	};

private:
	
	// Read Viewport dimension
	void getViewport() { glGetIntegerv( GL_VIEWPORT, vp ); };

};
typedef boost::shared_ptr<ImageBuffer> ImageBufferPtr;

class RecordingTools
{
public:
	RecordingTools();
	~RecordingTools();

	void init();
	void save_pixels();
	void writeTGA();
    void saveFrame();

	void setFilenamePrefix( std::string prefix ) { m_sPrefix = prefix; };

	void startRecording() { m_record = true;  };
	void stopRecording()  { m_record = false; };
	void toggleRecording() 
    { 
        m_record = !m_record;
        if( !m_record )
        { 
            m_buffer_ready = true; 
        }
    };
	bool isRecording()    { return m_record;  };

	int currentFrame, maxFrame;

	bool isWritingFiles() { return m_writeFiles; };
	int getWritingIndex() { return m_writeIndex; };

    bool buffer_ready() { return m_buffer_ready; };

private:

	std::queue<ImageBufferPtr> m_images;
	std::string m_sPrefix;

	bool m_record;
	bool m_writeFiles;
    bool m_buffer_ready;
	int  m_writeIndex;

};
typedef boost::shared_ptr<RecordingTools> RecordingToolsPtr;

#endif
