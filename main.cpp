

#include <iostream>
#include <vector>
#include <random>


#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "tools/stb_image.h"
#include "tools/stb_image_write.h"

#include "glm-master/glm/glm.hpp"


struct float4 {
	glm::vec4 real;
	glm::vec4 imag;
};

struct uchar4 {
	unsigned char x, y, z, w;
};

unsigned int reverseBits(unsigned int num, int log2n)
{
	unsigned int reversedNum = 0;
	for (int i = 0; i < log2n; ++i) {
		reversedNum <<= 1;
		reversedNum |= (num & 1);
		num >>= 1;
	}
	return reversedNum;
}


void fftx(std::vector<float4>& data, int width) {
	int N = width;
	int M = std::log2(N);

	for (int l = 0; l < width; l++) {
		for (unsigned int i = 0; i < N; i++) {
			unsigned int j = reverseBits(i, M);
			if (j > i) {
				std::swap(data[i + l * width], data[j + l * width]);
			}
		}
	}

	for (int l = 0; l < width; l++) {
		std::vector<float4> line(width);
		for (unsigned int i = 0; i < N; i++) {
			line[i] = data[i + l * width];
		}

		for (int k = 1; k <= M; k++) {
			int m = 1 << k;
			for (int j = 0; j < N; j += m) {
				for (int i = 0; i < m / 2; i++) {
					float theta = -2.0f * 3.141592f * (i - N / 2) / m;
					float fsin = std::sinf(theta);
					float fcos = std::cosf(theta);
					float4 t;
					float4 u;
					t.real = line[j + i + m / 2].real * fcos - line[j + i + m / 2].imag * fsin;
					t.imag = line[j + i + m / 2].real * fsin + line[j + i + m / 2].imag * fcos;
					u.real = line[j + i].real;
					u.imag = line[j + i].imag;
					line[j + i].real = u.real + t.real;
					line[j + i].imag = u.imag + t.imag;
					line[j + i + m / 2].real = u.real - t.real;
					line[j + i + m / 2].imag = u.imag - t.imag;
				}
			}
		}

		for (int n = 0; n < width; n++) {
			if (n % 2 == 0) {
				data[l * width + n].real = line[n].real;
				data[l * width + n].imag = line[n].imag;
			}
			else {
				data[l * width + n].real = line[n].real;
				data[l * width + n].imag = line[n].imag;
			}
		}
	}
}

void ffty(std::vector<float4>& data, int width) {
	int N = width;
	int M = std::log2(N);

	for (int l = 0; l < width; l++) {
		for (unsigned int i = 0; i < N; i++) {
			unsigned int j = reverseBits(i, M);
			if (j > i) {
				std::swap(data[(i * width) + l], data[j * width + l]);
			}
		}
	}

	for (int l = 0; l < width; l++) {
		std::vector<float4> line(width);
		for (unsigned int i = 0; i < N; i++) {
			line[i] = data[i * width + l];
		}


		for (int k = 1; k <= M; k++) {
			int m = 1 << k;
			for (int j = 0; j < N; j += m) {
				for (int i = 0; i < m / 2; i++) {
					float theta = -2.0f * 3.141592f * (i - N / 2) / m;
					float fsin = std::sinf(theta);
					float fcos = std::cosf(theta);
					float4 t;
					float4 u;
					t.real = line[j + i + m / 2].real * fcos - line[j + i + m / 2].imag * fsin;
					t.imag = line[j + i + m / 2].real * fsin + line[j + i + m / 2].imag * fcos;
					u.real = line[j + i].real;
					u.imag = line[j + i].imag;
					line[j + i].real = u.real + t.real;
					line[j + i].imag = u.imag + t.imag;
					line[j + i + m / 2].real = u.real - t.real;
					line[j + i + m / 2].imag = u.imag - t.imag;
				}
			}
		}

		for (int n = 0; n < width; n++) {
			if (n % 2 == 0) {
				data[l + width * n].real = line[n].real;
				data[l + width * n].imag = line[n].imag;
			}
			else {
				data[l + width * n].real = line[n].real;
				data[l + width * n].imag = line[n].imag;
			}
		}
	}
}

void ifftx(std::vector<float4>& data, int width) {
	int N = width;
	int M = std::log2(N);
	for (int l = 0; l < width; l++) {
		for (unsigned int i = 0; i < N; i++) {
			unsigned int j = reverseBits(i, M);
			if (j > i) {
				std::swap(data[i + width * l], data[j + width * l]);
			}
		}
	}

	for (int l = 0; l < width; l++) {
		std::vector<float4> line(width);
		for (unsigned int i = 0; i < N; i++) {
			line[i] = data[i + width * l];
		}


		for (int k = 1; k <= M; k++) {
			int m = 1 << k;
			for (int j = 0; j < N; j += m) {
				for (int i = 0; i < m / 2; i++) {
					float theta = 2.0f * 3.141592f * (i - N) / m;
					float fsin = std::sinf(theta);
					float fcos = std::cosf(theta);
					float4 t;
					float4 u;
					t.real = line[j + i + m / 2].real * fcos - line[j + i + m / 2].imag * fsin;
					t.imag = line[j + i + m / 2].real * fsin + line[j + i + m / 2].imag * fcos;
					u.real = line[j + i].real;
					u.imag = line[j + i].imag;
					line[j + i].real = u.real + t.real;
					line[j + i].imag = u.imag + t.imag;
					line[j + i + m / 2].real = u.real - t.real;
					line[j + i + m / 2].imag = u.imag - t.imag;
				}
			}
		}

		for (int n = 0; n < width; n++) {
			if (n % 2 == 0) {
				data[l * width + n].real = line[n].real / (float)N;
				data[l * width + n].imag = line[n].imag / (float)N;
			}
			else {
				data[l * width + n].real = -line[n].real / (float)N;
				data[l * width + n].imag = -line[n].imag / (float)N;
			}
		}
	}
}


void iffty(std::vector<float4>& data, int width) {
	int N = width;
	int M = std::log2(N);
	for (int l = 0; l < width; l++) {
		for (unsigned int i = 0; i < N; i++) {
			unsigned int j = reverseBits(i, M);
			if (j > i) {
				std::swap(data[i * width + l], data[j * width + l]);
			}
		}
	}

	for (int l = 0; l < width; l++) {
		std::vector<float4> line(width);
		for (unsigned int i = 0; i < N; i++) {
			line[i] = data[i * width + l];
		}


		for (int k = 1; k <= M; k++) {
			int m = 1 << k;
			for (int j = 0; j < N; j += m) {
				for (int i = 0; i < m / 2; i++) {
					float theta = 2.0f * 3.141592f * (i - N) / m;
					float fsin = std::sinf(theta);
					float fcos = std::cosf(theta);
					float4 t;
					float4 u;
					t.real = line[j + i + m / 2].real * fcos - line[j + i + m / 2].imag * fsin;
					t.imag = line[j + i + m / 2].real * fsin + line[j + i + m / 2].imag * fcos;
					u.real = line[j + i].real;
					u.imag = line[j + i].imag;
					line[j + i].real = u.real + t.real;
					line[j + i].imag = u.imag + t.imag;
					line[j + i + m / 2].real = u.real - t.real;
					line[j + i + m / 2].imag = u.imag - t.imag;
				}
			}
		}

		for (int n = 0; n < width; n++) {
			if (n % 2 == 0) {
				data[l + width * n].real = line[n].real / (float)N;
				data[l + width * n].imag = line[n].imag / (float)N;
			}
			else {
				data[l + width * n].real = -line[n].real / (float)N;
				data[l + width * n].imag = -line[n].imag / (float)N;
			}
		}
	}
}

int main()
{
	int width, height, channels;
	unsigned char* image = stbi_load("textures/test_out.png", &width, &height, &channels, 4);
	if (image == nullptr)
	{
		std::cout << "Error loading image" << std::endl;
		return 1;
	}
	std::cout << "Image loaded: " << width << "x" << height << "x" << channels << std::endl;

	std::random_device rnd;

	std::vector<float4> pixels(width * height);
	std::vector<float4> src(width * width);
	
	for (int i = 0; i < width * height; i++) {
		pixels[i].real.x = image[i * 4 + 0] / 255.0f;
		pixels[i].real.y = image[i * 4 + 1] / 255.0f;
		pixels[i].real.z = image[i * 4 + 2] / 255.0f;
		pixels[i].real.w = image[i * 4 + 3] / 255.0f;
		pixels[i].imag.x = 0.0f;
		pixels[i].imag.y = 0.0f;
		pixels[i].imag.z = 0.0f;
		pixels[i].imag.w = 0.0f;
		src[i] = pixels[i];
	}

	stbi_image_free(image);

	int blurWidth, blurHeight;
	image = stbi_load("textures/blur.png", &blurWidth, &blurHeight, &channels, 4);
	if (image == nullptr)
	{
		std::cout << "Error loading image" << std::endl;
		return 1;
	}
	std::cout << "Image loaded: " << blurWidth << "x" << blurHeight << "x" << channels << std::endl;


	std::vector<float4> blur(blurWidth * blurHeight);

	for (int i = 0; i < blurWidth * blurHeight; i++) {
		blur[i].real.x = image[i * 4 + 0] / 255.0f;
		blur[i].real.y = image[i * 4 + 1] / 255.0f;
		blur[i].real.z = image[i * 4 + 2] / 255.0f;
		blur[i].real.w = image[i * 4 + 3] / 255.0f;
		blur[i].imag.x = 0.0f;
		blur[i].imag.y = 0.0f;
		blur[i].imag.z = 0.0f;
		blur[i].imag.w = 0.0f;
	}

	stbi_image_free(image);

	fftx(blur, blurWidth);

	std::vector<uchar4> outPixels(width * height);
	for (int i = 0; i < blurWidth * blurHeight; i++) {
		outPixels[i].x = std::max(0.0f, std::min(blur[i].real.x * 1, 255.0f));
		outPixels[i].y = std::max(0.0f, std::min(blur[i].real.y * 1, 255.0f));
		outPixels[i].z = std::max(0.0f, std::min(blur[i].real.z * 1, 255.0f));
		outPixels[i].w = std::max(0.0f, std::min(blur[i].real.w * 1, 255.0f));
		outPixels[i].w = 255;
	}

	stbi_write_png("textures/test_blur_out0.png", blurWidth, blurHeight, 4, outPixels.data(), 0);

	ffty(blur, blurHeight);
	for (int i = 0; i < blurWidth * blurHeight; i++) {
		outPixels[i].x = std::max(0.0f, std::min(blur[i].real.x * 1, 255.0f));
		outPixels[i].y = std::max(0.0f, std::min(blur[i].real.y * 1, 255.0f));
		outPixels[i].z = std::max(0.0f, std::min(blur[i].real.z * 1, 255.0f));
		outPixels[i].w = std::max(0.0f, std::min(blur[i].real.w * 1, 255.0f));
		outPixels[i].w = 255;
	}

	stbi_write_png("textures/test_blur_out1.png", blurWidth, blurHeight, 4, outPixels.data(), 0);


	fftx(pixels, width);


	for (int i = 0; i < width * height; i++) {
		outPixels[i].x = std::max(0.0f, std::min(pixels[i].real.x * 1, 255.0f));
		outPixels[i].y = std::max(0.0f, std::min(pixels[i].real.y * 1, 255.0f));
		outPixels[i].z = std::max(0.0f, std::min(pixels[i].real.z * 1, 255.0f));
		outPixels[i].w = std::max(0.0f, std::min(pixels[i].real.w * 1, 255.0f));
		outPixels[i].w = 255;
	}

	stbi_write_png("textures/test_out0.png", width, height, 4, outPixels.data(), 0);

	ffty(pixels, width);

	for (int i = 0; i < width * height; i++) {
		outPixels[i].x = std::max(0.0f, std::min(pixels[i].real.x * 1, 255.0f));
		outPixels[i].y = std::max(0.0f, std::min(pixels[i].real.y * 1, 255.0f));
		outPixels[i].z = std::max(0.0f, std::min(pixels[i].real.z * 1, 255.0f));
		outPixels[i].w = std::max(0.0f, std::min(pixels[i].real.w * 1, 255.0f));
		outPixels[i].w = 255;
	}

	stbi_write_png("textures/test_out1.png", width, height, 4, outPixels.data(), 0);


	//dftx(blur, width);
	//dfty(blur, width);

	std::vector<float4> tempBuf(width * height);

	memcpy_s(tempBuf.data(), sizeof(float4) * width * height, pixels.data(), sizeof(float4) * width * height);

	for (int i = 0; i < width * height; i++) {
		float4 sum;
		memset(&sum, 0, sizeof(float4));
		sum.real.x += tempBuf[i].real.x * blur[i].real.x - tempBuf[i].imag.x * blur[i].imag.x;
		sum.real.y += tempBuf[i].real.y * blur[i].real.y - tempBuf[i].imag.y * blur[i].imag.y;
		sum.real.z += tempBuf[i].real.z * blur[i].real.z - tempBuf[i].imag.z * blur[i].imag.z;
		sum.real.w += tempBuf[i].real.w * blur[i].real.w - tempBuf[i].imag.w * blur[i].imag.w;
		sum.imag.x += tempBuf[i].real.x * blur[i].imag.x + tempBuf[i].imag.x * blur[i].real.x;
		sum.imag.y += tempBuf[i].real.y * blur[i].imag.y + tempBuf[i].imag.y * blur[i].real.y;
		sum.imag.z += tempBuf[i].real.z * blur[i].imag.z + tempBuf[i].imag.z * blur[i].real.z;
		sum.imag.w += tempBuf[i].real.w * blur[i].imag.w + tempBuf[i].imag.w * blur[i].real.w;

	}


	for (int i = 0; i < width * height; i++) {
		outPixels[i].x = std::max(0.0f, std::min(pixels[i].real.x * 255, 255.0f));
		outPixels[i].y = std::max(0.0f, std::min(pixels[i].real.y * 255, 255.0f));
		outPixels[i].z = std::max(0.0f, std::min(pixels[i].real.z * 255, 255.0f));
		outPixels[i].w = std::max(0.0f, std::min(pixels[i].real.w * 255, 255.0f));
	}

	stbi_write_png("textures/test_out2.png", width, height, 4, outPixels.data(), 0);

	iffty(pixels, width);

	for (int i = 0; i < width * height; i++) {
		outPixels[i].x = std::max(0.0f, std::min(pixels[i].real.x * 255, 255.0f));
		outPixels[i].y = std::max(0.0f, std::min(pixels[i].real.y * 255, 255.0f));
		outPixels[i].z = std::max(0.0f, std::min(pixels[i].real.z * 255, 255.0f));
		outPixels[i].w = std::max(0.0f, std::min(pixels[i].real.w * 255, 255.0f));
	}

	stbi_write_png("textures/test_out3.png", width, height, 4, outPixels.data(), 0);

	ifftx(pixels, width);

	for (int i = 0; i < width * height; i++) {
		outPixels[i].x = std::max(0.0f, std::min(pixels[i].real.r * 255, 255.0f));
		outPixels[i].y = std::max(0.0f, std::min(pixels[i].real.g * 255, 255.0f));
		outPixels[i].z = std::max(0.0f, std::min(pixels[i].real.b * 255, 255.0f));
		outPixels[i].w = std::max(0.0f, std::min(pixels[i].real.a * 255, 255.0f));
	}

	stbi_write_png("textures/test_out4.png", width, height, 4, outPixels.data(), 0);

	return 0;
} 