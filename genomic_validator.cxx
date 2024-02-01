#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

// Function prototypes
bool validateFormat(std::string const &fileName);
bool validateHeaderLine(std::string_view line);
bool checkDataLines(std::string_view line);

int main(int argc, char *argv[])
{
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <VCF filename>\n";
		return EXIT_FAILURE;
	}

	std::string fileName = argv[1];
	if (!validateFormat(fileName)) {
		std::cerr << "Invalid VCF file format.\n";
		return EXIT_FAILURE;
	}

	std::cout << "VCF file is valid.\n";
	return EXIT_SUCCESS;
}

bool isValidAlt(std::string_view alt)
{
	// Regular expression to match valid ALT alleles: bases, '*', or symbolic alleles (angle brackets not covered here)
	static std::regex const altRegex("^([ACGTN*]+|<[^>]+>)(,[ACGTN*]+|,<[^>]+>)*$");
	return std::regex_match(alt.cbegin(), alt.cend(), altRegex);
}

bool validateFormat(std::string const &fileName)
{
	std::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
	if (!file.is_open()) {
		std::cerr << "Failed to open file: " << fileName << '\n';
		return false;
	}

	boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
	if (!fileName.ends_with(".vcf")) {
		in.push(boost::iostreams::gzip_decompressor());
	}
	in.push(file);

	std::istream inf(&in);
	std::string line;
	bool headerLineFound = false;
	while (std::getline(inf, line)) {
		if (line.starts_with("##")) { // Meta-information lines
			if (!validateHeaderLine(line)) {
				return false;
			}
		} else if (line.starts_with("#") && !headerLineFound) { // Column header line
			headerLineFound = true;
			// Optional: Validate the content of the header line
			// if (!validateHeaderColumns(line)) return false;
		} else if (headerLineFound) { // Data lines
			// Process data lines
			if (!checkDataLines(line)) {
				return false;
			}
		} else {
			std::cerr << "Unexpected line format: " << line << '\n';
			return false;
		}
	}

	if (!headerLineFound) {
		std::cerr << "Missing column header line.\n";
		return false;
	}

	return true;
}

bool validateFileFormatLine(std::string_view line)
{
	static std::regex const fileFormatRegex("##fileformat=VCFv(\\d+\\.\\d+)");
	if (!std::regex_match(line.cbegin(), line.cend(), fileFormatRegex)) {
		std::cerr << "Invalid file format version: " << line << '\n';
		return false;
	}
	return true;
}

bool validateContigLine(std::string_view line)
{
	static std::regex const contigRegex("##contig=<ID=[^,]+(,length=\\d+)?(,.*)?>");
	if (!std::regex_match(line.cbegin(), line.cend(), contigRegex)) {
		std::cerr << "Invalid contig line: " << line << '\n';
		return false;
	}
	return true;
}

bool validateAltLine(std::string_view line)
{
	static std::regex const altRegex("##ALT=<ID=[^,]+,Description=\"[^\"]+\">");
	if (!std::regex_match(line.cbegin(), line.cend(), altRegex)) {
		std::cerr << "Invalid ALT line: " << line << '\n';
		return false;
	}
	return true;
}

bool validateSampleOrPedigreeLine(std::string_view line)
{
	// Basic structure check; can be more specific based on VCF version and use case
	if (line.starts_with("##SAMPLE=") || line.starts_with("##PEDIGREE=")) {
		return true; // Assuming well-formed for this example
	}
	std::cerr << "Invalid SAMPLE or PEDIGREE line: " << line << '\n';
	return false;
}

bool validateHeaderLine(std::string_view line)
{
	static std::regex const infoFormatRegex(
		"##(INFO|FORMAT)=<"
		"ID=[^,]+,"
		"Number=([\\.\\dAGRU]|-?\\d+),"
		"Type=(Integer|Float|Flag|Character|String),"
		"Description=\"[^\"]+\""
		"(,[^,]+=\"[^\"]+\")*>");
	static std::regex const filterRegex(
		"##FILTER=<"
		"ID=[^,]+,"
		"Description=\"[^\"]+\""
		">");

	if (line.starts_with("##INFO=") || line.starts_with("##FORMAT=")) {
		if (!std::regex_match(line.cbegin(), line.cend(), infoFormatRegex)) {
			std::cerr << "Invalid INFO or FORMAT line: " << line << '\n';
			return false;
		}
		return true;
	} else if (line.starts_with("##FILTER=")) {
		if (!std::regex_match(line.cbegin(), line.cend(), filterRegex)) {
			std::cerr << "Invalid FILTER line: " << line << '\n';
			return false;
		}
		return true;
	} else if (line.starts_with("##fileformat=")) {
		return validateFileFormatLine(line);
	} else if (line.starts_with("##contig=")) {
		return validateContigLine(line);
	} else if (line.starts_with("##ALT=")) {
		return validateAltLine(line);
	} else if (line.starts_with("##SAMPLE=") || line.starts_with("##PEDIGREE=")) {
		return validateSampleOrPedigreeLine(line);
	}

	if (line.starts_with("##")) {
		// Optionally, perform a basic structure check here, if necessary
		return true; // Accepts any well-formed header lines not covered by specific checks
	}

	std::cerr << "Unknown header format: " << line << '\n';
	return false;
}

bool checkTitleLine(std::string const &line)
{
	std::istringstream iss(line);
	std::vector<std::string> columns;
	std::string column;

	while (iss >> column) {
		columns.push_back(column);
	}

	// Check for required columns
	std::vector<std::string> const requiredColumns = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};
	if (columns.size() < requiredColumns.size()) {
		std::cerr << "Insufficient columns in title line.\n";
		return false;
	}

	return std::equal(columns.begin(), columns.begin() + requiredColumns.size(), requiredColumns.begin());
}

bool checkHeader(std::string const &line)
{
	if (line.starts_with("##fileformat=")) {
		return validateHeaderLine(line);
	} else if (line.starts_with("##INFO=") || line.starts_with("##FILTER=") || line.starts_with("##FORMAT=")) {
		// Add specific checks for each type of meta-information line if necessary
		return true;
	} else if (line[0] == '#' && line[1] != '#') { // Title line
		return checkTitleLine(line);
	} else {
		std::cerr << "Unknown header format: " << line << '\n';
		return false;
	}
}

std::vector<std::string_view> split(std::string_view str, char delimiter)
{
	std::vector<std::string_view> result;
	size_t start = 0;
	size_t end = 0;

	while ((end = str.find(delimiter, start)) != std::string_view::npos) {
		// Add a string_view of the token to the result
		result.emplace_back(str.substr(start, end - start));
		start = end + 1; // Move past the delimiter
	}
	// Add the last token after the final delimiter
	result.emplace_back(str.substr(start));

	return result;
}

bool isValidBase(std::string_view base)
{
	static std::regex const baseRegex("^[ACGTNacgtn]+$");
	return std::regex_match(base.cbegin(), base.cend(), baseRegex);
}

bool isValidGenotype(std::string_view gt)
{
	static std::regex const gtRegex("^(\\d+|\\.)([/|](\\d+|\\.))?$");
	return std::regex_match(gt.cbegin(), gt.cend(), gtRegex);
}

bool isNonNegativeInteger(std::string_view value)
{
	static std::regex const intRegex("^\\d+$");
	return std::regex_match(value.cbegin(), value.cend(), intRegex);
}

bool isListOfNonNegativeIntegers(std::string_view value)
{
	static std::regex const listRegex("^\\d+(,\\d+)*$");
	return std::regex_match(value.cbegin(), value.cend(), listRegex);
}

bool isFloat(std::string_view value)
{
	static std::regex const floatRegex("^[+-]?([0-9]*[.])?[0-9]+$");
	return std::regex_match(value.cbegin(), value.cend(), floatRegex);
}

bool isBoolean(std::string_view value)
{
	return value == "0" || value == "1";
}

bool isHumanChromosome(std::string_view chrom)
{
	// clang-format off
    static const std::set<std::string_view> humanChromosomesWithoutPrefix = {
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
        "11", "12", "13", "14", "15", "16", "17", "18", "19", 
        "20", "21", "22", "X", "Y", "MT"
    };

    static const std::set<std::string_view> humanChromosomesWithPrefix = {
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
        "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"
    };
	// clang-format on

	return humanChromosomesWithoutPrefix.find(chrom) != humanChromosomesWithoutPrefix.end()
		|| humanChromosomesWithPrefix.find(chrom) != humanChromosomesWithPrefix.end();
}

bool checkFormatAndSamples(std::vector<std::string_view> const &fields, size_t formatIndex)
{
	if (formatIndex >= fields.size()) {
		std::cerr << "FORMAT field missing or invalid\n";
		return false;
	}

	auto formatDescriptors = split(fields[formatIndex], ':');

	for (size_t i = formatIndex + 1; i < fields.size(); ++i) {
		auto sampleData = split(fields[i], ':');
		if (sampleData.size() != formatDescriptors.size()) {
			std::cerr << "Sample data does not match FORMAT descriptors\n";
			return false;
		}

		for (size_t j = 0; j < formatDescriptors.size(); ++j) {
			auto descriptor = formatDescriptors[j];
			auto data = sampleData[j];

			if (descriptor == "GT" && !isValidGenotype(data)) {
				std::cerr << "Invalid genotype data: " << data << '\n';
				return false;
			}

			if ((descriptor == "DP" || descriptor == "GQ") && !isNonNegativeInteger(data)) {
				std::cerr << "Invalid data for " << descriptor << ": " << data << '\n';
				return false;
			}

			if (descriptor == "AD" && !isListOfNonNegativeIntegers(data)) {
				std::cerr << "Invalid allele depth data: " << data << '\n';
				return false;
			}

			if (descriptor == "PL" && !isListOfNonNegativeIntegers(data)) {
				std::cerr << "Invalid phred-scaled genotype likelihoods data: " << data << '\n';
				return false;
			}

			if (descriptor == "MQ" && !isNonNegativeInteger(data)) {
				std::cerr << "Invalid mapping quality data: " << data << '\n';
				return false;
			}

			if (descriptor == "SB" && !isListOfNonNegativeIntegers(data)) {
				std::cerr << "Invalid strand bias data: " << data << '\n';
				return false;
			}

			if (descriptor == "MQ0" && !isNonNegativeInteger(data)) {
				std::cerr << "Invalid MQ0 data: " << data << '\n';
				return false;
			}

			if (descriptor == "HRun" && !isNonNegativeInteger(data)) {
				std::cerr << "Invalid homopolymer run length data: " << data << '\n';
				return false;
			}

			if (descriptor == "AF" && !isFloat(data)) {
				std::cerr << "Invalid allele frequency data: " << data << '\n';
				return false;
			}

			if (descriptor == "AC" && !isNonNegativeInteger(data)) {
				std::cerr << "Invalid allele count data: " << data << '\n';
				return false;
			}

			if (descriptor == "AN" && !isNonNegativeInteger(data)) {
				std::cerr << "Invalid total number of alleles data: " << data << '\n';
				return false;
			}

			if (descriptor == "BaseQRankSum" && !isFloat(data)) {
				std::cerr << "Invalid Base Quality Rank Sum Test data: " << data << '\n';
				return false;
			}

			if (descriptor == "ReadPosRankSum" && !isFloat(data)) {
				std::cerr << "Invalid Read Position Rank Sum Test data: " << data << '\n';
				return false;
			}

			if (descriptor == "FS" && !isFloat(data)) {
				std::cerr << "Invalid Fisher Strand Bias data: " << data << '\n';
				return false;
			}

			if (descriptor == "SOR" && !isFloat(data)) {
				std::cerr << "Invalid Strand Odds Ratio data: " << data << '\n';
				return false;
			}

			if (descriptor == "MQRankSum" && !isFloat(data)) {
				std::cerr << "Invalid Mapping Quality Rank Sum Test data: " << data << '\n';
				return false;
			}

			if (descriptor == "QD" && !isFloat(data)) {
				std::cerr << "Invalid Quality by Depth data: " << data << '\n';
				return false;
			}

			if (descriptor == "RPA" && !isListOfNonNegativeIntegers(data)) {
				std::cerr << "Invalid Repeat unit number data: " << data << '\n';
				return false;
			}

			if (descriptor == "RU" && data.empty()) {
				std::cerr << "Invalid Repeat unit data: " << data << '\n';
				return false;
			}

			if (descriptor == "STR" && !isBoolean(data)) {
				std::cerr << "Invalid Short Tandem Repeat data: " << data << '\n';
				return false;
			}

			// Additional checks for other descriptors can be added here
		}
	}

	return true;
}

int stringViewToInt(std::string_view sv)
{
	int value = 0;
	auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), value);
	if (ec == std::errc::invalid_argument || ec == std::errc::result_out_of_range) {
		throw std::runtime_error("Conversion error");
	}
	return value;
}

float stringViewToFloat(std::string_view sv)
{
	float value = 0;
	auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), value);
	if (ec == std::errc::invalid_argument || ec == std::errc::result_out_of_range) {
		throw std::runtime_error("Conversion error");
	}
	return value;
}

bool checkDataLines(std::string_view line)
{
	auto fields = split(line, '\t');

	// Basic check for the number of fields
	size_t const expectedFieldCount = 8;
	if (fields.size() < expectedFieldCount) {
		std::cerr << "Invalid data line (not enough fields): " << line << '\n';
		return false;
	}

	// Validate CHROM - simple check for non-empty string
	if (fields[0].empty()) {
		std::cerr << "Invalid CHROM field: " << fields[0] << '\n';
		return false;
	}

	// Check if CHROM field is a human chromosome
	if (!isHumanChromosome(fields[0])) {
		std::cerr << "Non-human chromosome found: " << fields[0] << '\n';
		return false;
	}

	// Validate POS - should be a positive integer
	try {
		int pos = stringViewToInt(fields[1]);
		if (pos <= 0) {
			std::cerr << "Invalid POS field: " << fields[1] << '\n';
			return false;
		}
	} catch (std::invalid_argument &e) {
		std::cerr << "Invalid POS field (not an integer): " << fields[1] << '\n';
		return false;
	}

	// Validate ID - should be a string or '.'
	if (fields[2] != "." && fields[2].empty()) {
		std::cerr << "Invalid ID field: " << fields[2] << '\n';
		return false;
	}

	// Validate REF - should be one of A, C, G, T, N
	if (!isValidBase(fields[3])) {
		std::cerr << "Invalid REF field: " << fields[3] << '\n';
		return false;
	}

	// Validate ALT field
	if (!isValidAlt(fields[4])) { // Assuming ALT is the fifth column (0-based indexing)
		std::cerr << "Invalid ALT field: " << fields[4] << '\n';
		return false;
	}

	// Validate QUAL - should be a float or '.'
	if (fields[5] != ".") {
		try {
			float qual = stringViewToFloat(fields[5]);
			if (qual < 0) {
				std::cerr << "Invalid QUAL field: " << fields[5] << '\n';
				return false;
			}
		} catch (std::invalid_argument &e) {
			std::cerr << "Invalid QUAL field (not a float): " << fields[5] << '\n';
			return false;
		}
	}

	// Validate FILTER - should be a string or '.'
	if (fields[6] != "." && fields[6].empty()) {
		std::cerr << "Invalid FILTER field: " << fields[6] << '\n';
		return false;
	}

	// Validate INFO - additional information in key=value format; complex validation can be added here
	if (fields[7].empty()) {
		std::cerr << "Invalid INFO field: " << fields[7] << '\n';
		return false;
	}

	// Check FORMAT and sample-specific columns
	size_t const formatFieldIndex = 8; // FORMAT field is the 9th column (0-based index)
	if (!checkFormatAndSamples(fields, formatFieldIndex)) {
		return false;
	}

	return true;
}
