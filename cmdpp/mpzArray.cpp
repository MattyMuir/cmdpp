#include "mpzArray.h"

mpzArray::mpzArray()
{
	mData = std::vector<std::shared_ptr<mpzxx>>();
	size = 0;
}

mpzArray::mpzArray(int count)
{
	mData = std::vector<std::shared_ptr<mpzxx>>();
	size = 0;

	for (int i = 0; i < count; i++)
		this->Back();
}

mpzArray::~mpzArray()
{
	this->Clear();
}

void mpzArray::Reserve(int count)
{
	int originalSize = size;
	if (count < mData.size())
	{
		originalSize = 0;
		this->Clear();
	}
	while (mData.size() < count)
		this->Back();
	size = originalSize;
}

void mpzArray::Clear()
{
	mData.clear();
}

void mpzArray::TempClear()
{
	size = 0;

#ifdef _DEBUG
	for (int i = 0; i < mData.size(); i++)
	{
		mpz_set_ui((*mData[i]).data, 0);
	}
#endif
}

mpz_t& mpzArray::operator[](int i)
{
#ifdef _DEBUG
	if (i >= size) { throw; }
#endif

	return (*mData[i]).data;
}

int mpzArray::Size()
{
	return size;
}

void mpzArray::Back()
{
	if (mData.size() <= size)
	{
		mData.push_back(std::make_shared<mpzxx>());
		mpz_init((*mData[mData.size() - 1]).data);
	}
	size++;
}

void mpzArray::Back(mpz_t element)
{
	if (mData.size() <= size)
	{
		mData.push_back(std::make_shared<mpzxx>());
		mpz_init_set((*mData[mData.size() - 1]).data, element);
	}
	else
	{
		mpz_set((*mData[size]).data, element);
	}
	size++;
}

void mpzArray::BackUI(uint64_t element)
{
	if (mData.size() <= size)
	{
		mData.push_back(std::make_shared<mpzxx>());
		mpz_init_set_ui((*mData[mData.size() - 1]).data, element);
	}
	else
	{
		mpz_set_ui((*mData[size]).data, element);
	}
	size++;
}

void mpzArray::DeepCopy(mpzArray& other)
{
	this->TempClear();
	for (int i = 0; i < other.Size(); i++)
	{
		this->Back(other[i]);
	}
}

void mpzArray::Resize(int newSize)
{
	mData.resize(newSize);
}