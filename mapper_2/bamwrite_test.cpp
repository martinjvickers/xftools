#include "common.h"

int main(int argc, char const ** argv)
{

    // Initialize BamFile, build header.
    seqan::BamFileOut bamIO(toCString("test_meh.bam"));

    seqan::BamHeader header;
    assignValueById(contigLengths(context(bamIO)), nameToId(contigNamesCache(context(bamIO)), "REFERENCE"), 10000);
    resize(header, 2);
    resize(header[0].tags, 2);
    header[0].type = seqan::BAM_HEADER_FIRST;
    header[0].tags[0].i1 = "VN";
    header[0].tags[0].i2 = "1.3";
    header[0].tags[1].i1 = "SO";
    header[0].tags[1].i2 = "coordinate";
    resize(header[1].tags, 2);
    header[1].type = seqan::BAM_HEADER_REFERENCE;
    header[1].tags[0].i1 = "SN";
    header[1].tags[0].i2 = "REFERENCE";
    header[1].tags[1].i1 = "LN";
    header[1].tags[1].i2 = "10000";
    writeHeader(bamIO, header);

    // Construct first records.
    seqan::BamAlignmentRecord record;

    record.qName = "READ0";
    record.flag = 2;
    record.rID = 0;
    record.beginPos = 0;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    writeRecord(bamIO, record);

    record.qName = "READ0";
    record.flag = 1;
    record.rID = 0;
    record.beginPos = 1;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    writeRecord(bamIO, record);

    record.qName = "READ0";
    record.flag = 3;
    record.rID = 0;
    record.beginPos = 2;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = seqan::BamAlignmentRecord::INVALID_REFID;
    record.pNext = seqan::BamAlignmentRecord::INVALID_POS;
    record.tLen = seqan::BamAlignmentRecord::INVALID_LEN;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    writeRecord(bamIO, record);

    // Force writing of everything.
    close(bamIO);
}
