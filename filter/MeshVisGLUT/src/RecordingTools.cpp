#include "RecordingTools.h"

RecordingTools::RecordingTools() : 
	m_record( false ), 
	m_writeFiles(false),
    m_buffer_ready(false),
	m_writeIndex(0),
	m_sPrefix("")
{
	currentFrame = 0;
	maxFrame = MAX_FRAMES_BEFORE_FLUSH;
};

RecordingTools::~RecordingTools()
{
};

void RecordingTools::init()
{
};

// read framebuffer and store im mem
void RecordingTools::save_pixels( )
{
	/* 	GL_PACK_ALIGNMENT:
		Specifies the alignment requirements for the start of each pixel row in memory. 
		The allowable values are
		1 (byte-alignment),
		2 (rows aligned to even-numbered bytes),
		4 (word-alignment), and
		8 (rows start on double-word boundaries). */

	glPixelStorei(GL_PACK_ALIGNMENT,1);

	// The screenshot data
	ImageBufferPtr buffer( new ImageBuffer() );
   
	// alloc mem
	buffer->init();

	// read pixel buffer
	buffer->readPixel();

	// store
	m_images.push(buffer);

	currentFrame++;

    m_buffer_ready = (currentFrame > maxFrame);
};

void RecordingTools::writeTGA()
{
	m_writeFiles = true;
	//m_writeIndex = 0;

	//maxFrame += MAX_FRAMES_BEFORE_FLUSH;
	
	FILE *ftga;
    char cFileName[64];
   
    // write buffered image-data to files
	while( !m_images.empty() )
	{
        ImageBufferPtr img = m_images.front();

		// generate the tga-header on-the-fly
		unsigned char header[6] = {
			img->width & 0x00FF,
			( img->width & 0xFF00) / 256.0,
			img->height & 0x00FF,
			(img->height & 0xFF00) / 256.0,
			24, 0	
		};

        bool stop = false;
        do 
		{
           sprintf(cFileName,"C:\\Frames\\%04d.tga", m_writeIndex++);
           ftga = fopen( cFileName, "rb" );
           stop = (ftga == NULL);
           fclose( ftga );
		}
		while( !stop );

		ftga = fopen(cFileName, "wb");
    
		fwrite(TGAheader, sizeof(unsigned char), 12, ftga);
		fwrite(header, sizeof(unsigned char), 6, ftga);
		fwrite(img->data, sizeof(GLubyte), img->image_size, ftga);
		
		fclose(ftga);

		m_images.pop();
	}

	m_writeFiles = false;
	//m_buffer_ready = false;
};

void RecordingTools::saveFrame()
{
	glPixelStorei(GL_PACK_ALIGNMENT,1);

	// Create buffer and read pixel
	ImageBufferPtr buffer( new ImageBuffer() );

    buffer->init();
	buffer->readPixel();

    // write Frame
    FILE *ftga;
    char cFileName[64];
   
	// generate the tga-header on-the-fly
	unsigned char header[6] = {
		buffer->width & 0x00FF,
		( buffer->width & 0xFF00) / 256.0,
		buffer->height & 0x00FF,
		( buffer->height & 0xFF00) / 256.0,
		24, 0	
	};

  do
	{
#ifdef _WIN32
    sprintf(cFileName,"C:\\Frames\\%04d.tga", m_writeIndex++);
#else
    sprintf(cFileName,"./Frames/%04d.tga", m_writeIndex++);
#endif
	}
	while( fopen( cFileName, "rb" ) != NULL);

	ftga = fopen(cFileName, "wb");
    
	fwrite(TGAheader, sizeof(unsigned char), 12, ftga);
	fwrite(header, sizeof(unsigned char), 6, ftga);
	fwrite(buffer->data, sizeof(GLubyte), buffer->image_size, ftga);
		
	fclose(ftga);

};
