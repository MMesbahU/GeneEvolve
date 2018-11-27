#include "format_vcf.h"


// vcf_st.data = nhap*nsnp, with 0=REF=false, 1=true=ALT. In C, 0 means false
bool format_vcf::write_vcf_file(std::string file_out_name, vcf_structure &vcf_st)
{
    std::ofstream outfile;
    outfile.open(file_out_name.c_str());

    if(!outfile)
    {
        std::cout << "Error: can not open the file [" + file_out_name + "] to read." << std::endl;
        return false;
    }

    unsigned nhap = vcf_st.data.size();
    unsigned nsnp = vcf_st.data[0].size();
    unsigned nind = vcf_st.SAMPLES.size();
    if (vcf_st.POS.size()!=nsnp)
    {
        std::cout << "Error: file [" + file_out_name + "] has problem." << std::endl;
        return false;
    }
    if (vcf_st.SAMPLES.size()*2!=nhap)
    {
        std::cout << "Error: file [" + file_out_name + "] has problem.2." << std::endl;
        return false;
    }

    // write meta lines, starts with ##
    for (unsigned i=0; i<vcf_st.meta_lines.size(); i++)
    {
        outfile << vcf_st.meta_lines[i] << std::endl;
    }

    // create head line, starts with #
    outfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (unsigned iind=0; iind<nind; iind++)
    {
        outfile << "\t" << vcf_st.SAMPLES[iind];
    }
    outfile << std::endl;


    // data
    for (unsigned isnp=0; isnp<nsnp; isnp++)
    {
        std::string str_QUAL=".";
        if (vcf_st.QUAL[isnp]!=-1)
        {
            str_QUAL = std::to_string(vcf_st.QUAL[isnp]);
        }
        outfile << vcf_st.CHROM[isnp] << "\t" << vcf_st.POS[isnp] << "\t" << vcf_st.ID[isnp] << "\t" <<
        vcf_st.REF[isnp] << "\t" << vcf_st.ALT[isnp] << "\t" << str_QUAL << "\t" << vcf_st.FILTER[isnp] << "\t" << vcf_st.INFO[isnp] << "\t" << vcf_st.FORMAT[isnp];
        for (int iind=0; iind<nind; iind++)
        {
            outfile << "\t" << vcf_st.data[2*iind][isnp] << "|" << vcf_st.data[2*iind+1][isnp];
        }
        outfile << std::endl;
    }

    outfile.close();

    return true;

}





// reads just biallelic SNPs
// vcf_st.data = nhap*nsnp, with 0=REF=false, 1=true=ALT. In C, 0 means false
bool format_vcf::read_vcf_file(std::string filename, vcf_structure &vcf_st)
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    // primary analysis

    inFile.setSiteOnly(true); // true, if you dont want to read genotypes
    vector<int> importIndexList;

    int numReadRecords = 0, numActualRecords=0;
    int failFilter=0, notBiallelic=0, inconsistent=0, otherChrom=0;

    if (!inFile.open(filename.c_str(), header))
    {
        std::cout << "\n Program could NOT open file : " << filename << std::endl;
        return false;
    }

    std::cout << "\n Loading VCF File     : " << filename << std::endl<<std::endl;
    int numSamples = header.getNumSamples();
    while (inFile.readRecord(record))
    {
        int flag=0;
        std::string cno       = record.getChromStr();
        int bp                = record.get1BasedPosition();
        std::string id        = record.getIDStr();
        std::string refAllele = record.getRefStr();
        std::string altAllele = record.getAltStr();
        float qual            = record.getQual();
        std::string filter    = record.getFilter().getString();
        std::string info      = ""; //record.getInfo();
        std::string currID    = cno + ":" + std::to_string(bp);
        if(id==".")
            id=currID;

        if (record.getNumAlts()>1) // nonBiallelic
        {
            notBiallelic++;
            flag = 1;
        }
        if (record.getFilter().getString(0).compare("PASS") != 0)
        {
            failFilter++;
        }
        char Rallele=0,Aallele=0;

        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
        {
            switch (refAllele[0])
            {
                case 'A': case 'a': Rallele = 1; break;
                case 'C': case 'c': Rallele = 2; break;
                case 'G': case 'g': Rallele = 3; break;
                case 'T': case 't': Rallele = 4; break;
                case 'D': case 'd': Rallele = 5; break;
                case 'I': case 'i': Rallele = 6; break;
                case 'R': case 'r': Rallele = 7; break;
                default:
                {
                    std::cout << " WARNING !!! Reference allele for SNP for [" << currID << "] is [" <<refAllele << "]. Will be ignored ..." << std::endl;
                    flag=1;
                    inconsistent++;
                }
            }
            if(flag==0)
            switch (altAllele[0])
            {
                case '0':  Aallele = 0; break;
                case 'A': case 'a': Aallele = 1; break;
                case 'C': case 'c': Aallele = 2; break;
                case 'G': case 'g': Aallele = 3; break;
                case 'T': case 't': Aallele = 4; break;
                case 'D': case 'd': Aallele = 5; break;
                case 'I': case 'i': Aallele = 6; break;
                case 'R': case 'r': Aallele = 7; break;
                default:
                {
                    std::cout << " WARNING !!! Alternate allele for SNP for [" << currID << "] is [" << altAllele << "]. Will be ignored ..." << std::endl;
                    flag=1;
                    inconsistent++;
                }
            }
        }
        else
        {
            Rallele = 7;
            if(strlen(refAllele.c_str())<strlen(altAllele.c_str()))
            Aallele=6;
            else
            Aallele=5;
        }


        if(flag==0)
        {
            ++numReadRecords;
        }

        numActualRecords++;
    }
    inFile.close();
    std::cout << " Number of true Markers read from VCF File                : " << numReadRecords << std::endl;

    vcf_st.data.resize(2*numSamples, std::vector<bool> (numReadRecords,false)); // nhap*nsnp
    vcf_st.CHROM.resize(numReadRecords);
    vcf_st.POS.resize(numReadRecords);
    vcf_st.ID.resize(numReadRecords);
    vcf_st.REF.resize(numReadRecords);
    vcf_st.ALT.resize(numReadRecords);
    vcf_st.QUAL.resize(numReadRecords);
    vcf_st.FILTER.resize(numReadRecords);
    vcf_st.INFO.resize(numReadRecords);
    vcf_st.FORMAT.resize(numReadRecords);


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    // secondary analysis

    //inFile.setSiteOnly(true); // if you dont want to read genotypes
    inFile.setSiteOnly(false);  // if you want to read genotypes
    numReadRecords = 0, numActualRecords=0;
    failFilter=0, notBiallelic=0, inconsistent=0, otherChrom=0;

    if (!inFile.open(filename.c_str(), header))
    {
        std::cout << "\n Program could NOT open file : " << filename << std::endl;
        return false;
    }


    numSamples = header.getNumSamples();
    vcf_st.SAMPLES.resize(numSamples);
    std::cout << " numSamples=" << numSamples << std::endl;
    for (int i=0; i<numSamples; i++)
    {
        vcf_st.SAMPLES[i] = header.getSampleName(i);
    }


    while (inFile.readRecord(record))
    {

        int flag=0;
        std::string cno       = record.getChromStr();
        int bp                = record.get1BasedPosition();
        std::string id        = record.getIDStr();
        std::string refAllele = record.getRefStr();
        std::string altAllele = record.getAltStr();
        float qual            = record.getQual();
        std::string filter    = record.getFilter().getString();
        std::string info      = ""; //record.getInfo();
        std::string format    = ""; //record.getFilter();
        std::string currID    = cno + ":" + std::to_string(bp);
        if(id==".")
            id=currID;


        if (record.getNumAlts()>1)
        {
            notBiallelic++;
            flag = 1;
        }
        if (record.getFilter().getString(0).compare("PASS") != 0)
        {
            failFilter++;
        }

        char Rallele=0,Aallele=0;
        if (strlen(refAllele.c_str()) == 1 && strlen(altAllele.c_str()) == 1)
        {
            switch (refAllele[0])
            {
                case 'A': case 'a': Rallele = 1; break;
                case 'C': case 'c': Rallele = 2; break;
                case 'G': case 'g': Rallele = 3; break;
                case 'T': case 't': Rallele = 4; break;
                case 'D': case 'd': Rallele = 5; break;
                case 'I': case 'i': Rallele = 6; break;
                case 'R': case 'r': Rallele = 7; break;
                default:
                {
                    std::cout << " WARNING !!! Reference allele for SNP for [" << currID << "] is [" <<refAllele << "]. Will be ignored ..." << std::endl;
                    flag=1;
                    inconsistent++;
                }
            }
            if(flag==0)
            switch (altAllele[0])
            {
                case '0':  Aallele = 0; break;
                case 'A': case 'a': Aallele = 1; break;
                case 'C': case 'c': Aallele = 2; break;
                case 'G': case 'g': Aallele = 3; break;
                case 'T': case 't': Aallele = 4; break;
                case 'D': case 'd': Aallele = 5; break;
                case 'I': case 'i': Aallele = 6; break;
                case 'R': case 'r': Aallele = 7; break;
                default:
                {
                    std::cout << " WARNING !!! Alternate allele for SNP for [" << currID << "] is [" << altAllele <<"]. Will be ignored ..." << std::endl;
                    flag=1;
                    inconsistent++;
                }
            }
        }
        else
        {
            Rallele = 7;
            if(strlen(refAllele.c_str()) < strlen(altAllele.c_str()))
            Aallele=6;
            else
            Aallele=5;
        }


        if(flag==0)
        {
            vcf_st.CHROM[numReadRecords]=cno;
            vcf_st.POS[numReadRecords]=bp;
            vcf_st.ID[numReadRecords]=id;
            vcf_st.REF[numReadRecords]=refAllele;
            vcf_st.ALT[numReadRecords]=altAllele;
            vcf_st.QUAL[numReadRecords]=qual;
            vcf_st.FILTER[numReadRecords]=filter;
            vcf_st.INFO[numReadRecords]=info;
            vcf_st.FORMAT[numReadRecords]=format;

            ////////////////////////////////////////////////
            // for samples
            for (int i = 0; i<numSamples; i++)
            {
                if(record.getNumGTs(i)==0)
                {
                    std::cout << "\n Empty Value for Individual : " << vcf_st.SAMPLES[i] << " at Marker : " <<
                    id << std::endl;
                    std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << std::endl;
                    //return false;
                }
                //if(record.getNumGTs(i)==1 && record.getChromStr()!="X")
                //{
                //    std::cout << "\n Single Autosomal Haplotype for Individual : " << individualName[i] << " at Marker : "
                //    << thisVariant.name << std::endl;
                //    std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << std::endl;
                //    return false;
                //}


                for (int j = 0; j<record.getNumGTs(i); j++) // for j=0,1. because it is biallelic
                {

                    int alleleIndex = record.getGT(i, j); // should be 0 or 1
                    //std::cout << "     record.getGT(i,j)=" << alleleIndex << std::endl;
                    vcf_st.data[2*i+j][numReadRecords] = (bool)alleleIndex; // data = nhap*nsnp with 0=REF=false, 1=true=ALT. In C, 0 means false
                }

            } // for samples

            ++numReadRecords;
        } // if (flag==0)
        numActualRecords++;
        importIndexList.push_back(flag);
    }// while

    std::cout << std::endl;
    std::cout << " Number of Markers read from VCF File                : " << numActualRecords << std::endl;
    std::cout << " Number of Markers with more than One Allele         : " << notBiallelic << std::endl;
    std::cout << " Number of Markers failing FILTER = PASS             : " << failFilter << std::endl;
    std::cout << " Number of Markers with inconsistent Ref/Alt Allele  : " << inconsistent << std::endl;
    std::cout << " Number of Markers on other chromosomes (Non-Ref)    : " << otherChrom << std::endl;

    if(numActualRecords==0)
    {
        std::cout << std::endl;
        std::cout << " No Markers recorded from VCF Input File : " << filename << std::endl;
        std::cout << " Please check the file properly.." << std::endl;
        std::cout << " Program Aborting ... " << std::endl;
        return false;
    }
    inFile.close();


    return true;

}


// output: std::vector<std::string> sample;
bool format_vcf::read_vcf_header_sample(std::string filename, std::vector<std::string> &sample)
{
    VcfFileReader inFile;
    VcfHeader header;

    inFile.setSiteOnly(true); // true, if you dont want to read genotypes


    if (!inFile.open(filename.c_str(), header))
    {
        std::cout << "\n Program could NOT open file : " << filename << std::endl;
        return false;
    }

    int numSamples = header.getNumSamples();

    sample.resize(numSamples);
    for (int i=0; i<numSamples; i++)
    {
        sample[i] = header.getSampleName(i);
    }

}
