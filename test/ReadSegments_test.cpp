/**
 * @file test/ReadSegments_test.cpp
 *
 * @date 2015-06-04
 *
 * @author Youri Hoogstrate
 *
 * @section LICENSE
 * segmentation-fold can predict RNA 2D structures including K-turns.
 * Copyright (C) 2012-2015 Youri Hoogstrate
 *
 * This file is part of segmentation-fold.
 *
 * segmentation-fold is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * segmentation-fold is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */



#define BOOST_TEST_MODULE ReadSegments



#include <iostream>

#include "config.hpp"

#include "Nucleotide.hpp"
#include "Pair.hpp"
#include "SubSequence.hpp"
#include "Sequence.hpp"

#include "Direction.hpp"
#include "Segment.hpp"
#include "SegmentTreeElement.hpp"
#include "SegmentTree.hpp"

#include "ReadSegments.hpp"



#include <boost/test/included/unit_test.hpp>

// See what happens if you remove this:
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>



BOOST_AUTO_TEST_SUITE(Testing_Segments)


/**
 * @brief checks whether parsing a custom file succeeds
 * 
 * @date 2015-06-04
 */
BOOST_AUTO_TEST_CASE(Test1)
{
	std::string filename = "tmp.readsegments_test_test1";
	
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile <<   "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
		   "<root>\n"
		   "	<motifs>\n"
		   "	</motifs>\n"
		   "</root>\n";
	myfile.close();
	
	
	SegmentTree segments;
	ReadSegments r = ReadSegments(filename, segments);
	
	BOOST_CHECK_EQUAL(segments.empty() , true);
	BOOST_CHECK_EQUAL(segments.size() , 0);
	
	unlink(filename.c_str());
}

/**
 * @brief checks whether parsing a custom file succeeds - testing with 5' -> 3'
 * 
 * @test
 * 
 * @date 2015-06-04
 */
BOOST_AUTO_TEST_CASE(Test2)
{
	std::string filename = "tmp.readsegments_test_test2";
	
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile <<   "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
		   "<root>\n"
		   "	<motifs>\n"
		   "		<motif><id>Kt-7 G2nA SAM riboswitch (H. marismortui)</id>\n"
		   "			<sequence_5prime>GAAGAA</sequence_5prime>\n"
		   "			<bonds          >   :::</bonds>\n"
		   "			<sequence_3prime>   AGG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions>\n"
		   "				<five_prime>true</five_prime>\n"
		   "				<three_prime>false</three_prime>\n"
		   "			</directions>\n"
		   "		</motif>\n"
		   "	</motifs>\n"
		   "</root>\n";
	myfile.close();
	
	SegmentTree segments;
	ReadSegments r = ReadSegments(filename, segments);
	
	BOOST_CHECK_EQUAL(segments.empty() , false);
	BOOST_CHECK_EQUAL(segments.size() , 1);
	
	Sequence sequence_5p = Sequence("GAAGAA");
	Sequence sequence_3p = Sequence("GGA");// Rotated; 5'->3'
	Position position_5p_s = sequence_5p.data.begin();
	Position position_5p_e = sequence_5p.data.end() - 1;
	Position position_3p_s = sequence_3p.data.begin();
	Position position_3p_e = sequence_3p.data.end() - 1;
	SubSequence subsequence_5p = SubSequence(position_5p_s, position_5p_e);
	SubSequence subsequence_3p = SubSequence(position_3p_s, position_3p_e);
	Segment *segment = segments.search(subsequence_5p, subsequence_3p);
	
	BOOST_REQUIRE(segment != NULL);
	BOOST_CHECK_EQUAL(segment->gibbs_free_energy , (float) - 11.1072);
	
	// Check bonds
	BOOST_CHECK_EQUAL(segment->bonds[0].first , 3);
	BOOST_CHECK_EQUAL(segment->bonds[0].second , 2);
	
	BOOST_CHECK_EQUAL(segment->bonds[1].first , 4);
	BOOST_CHECK_EQUAL(segment->bonds[1].second , 1);
	
	BOOST_CHECK_EQUAL(segment->bonds[2].first , 5);
	BOOST_CHECK_EQUAL(segment->bonds[2].second , 0);
	
	signed int i, j;
	bool s;
	
	for(unsigned int k = 0; k < 12; k++)								// Check tracing back bonds:
	{
		switch(k % 4)
		{
			case 0:
				BOOST_CHECK_EQUAL(segment->pop(i, j), true);
				BOOST_CHECK_EQUAL(i, -1);
				BOOST_CHECK_EQUAL(j, -1);
				break;
			case 1:
				BOOST_CHECK_EQUAL(segment->pop(i, j), true);
				BOOST_CHECK_EQUAL(i, -2);
				BOOST_CHECK_EQUAL(j, -2);
				break;
			case 2:
				BOOST_CHECK_EQUAL(segment->pop(i, j), true);
				BOOST_CHECK_EQUAL(i, -3);
				BOOST_CHECK_EQUAL(j, -3);
				break;
			case 3:
				BOOST_CHECK_EQUAL(segment->pop(i, j), false);
				break;
		}
	}
	
	unlink(filename.c_str());
}

/**
 * @brief checks whether parsing a custom file succeeds - testing with 3' -> 5'
 *
 * 5') GAAGAA (3'
 *        :::
 * 3')    AGG (5'
 *
 * rotates to
 *
 * 5') GGA    (3'
 *     :::
 * 3') AAGAAG (5'
 * 
 * @test
 * 
 * @date 2015-06-04
 * 
 * @todo also check segments pop() function
 */
BOOST_AUTO_TEST_CASE(Test3)
{
	std::string filename = "tmp.readsegments_test_test3";
	
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
		   "<root>\n"
		   "	<motifs>\n"
		   "		<motif><id>Kt-7 G2nA SAM riboswitch (H. marismortui)</id>\n"
		   "			<sequence_5prime>GAAGAA</sequence_5prime>\n"
		   "			<bonds          >   :::</bonds>\n"
		   "			<sequence_3prime>   AGG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions>\n"
		   "				<five_prime>false</five_prime>\n"
		   "				<three_prime>true</three_prime>\n"
		   "			</directions>\n"
		   "		</motif>\n"
		   "	</motifs>\n"
		   "</root>\n";
	myfile.close();
	
	SegmentTree segments;
	ReadSegments r = ReadSegments(filename, segments);
	
	
	BOOST_CHECK_EQUAL(segments.empty() , false);
	BOOST_CHECK_EQUAL(segments.size() , 1);
	
	Sequence sequence_5p = Sequence("GGA");
	Sequence sequence_3p = Sequence("GAAGAA");// Rotated; 5'->3'
	Position position_5p_s = sequence_5p.data.begin();
	Position position_5p_e = sequence_5p.data.end() - 1;
	Position position_3p_s = sequence_3p.data.begin();
	Position position_3p_e = sequence_3p.data.end() - 1;
	SubSequence subsequence_5p = SubSequence(position_5p_s, position_5p_e);
	SubSequence subsequence_3p = SubSequence(position_3p_s, position_3p_e);
	Segment *segment = segments.search(subsequence_5p, subsequence_3p);
	
	
	BOOST_REQUIRE(segment != NULL);
	BOOST_CHECK_EQUAL(segment->gibbs_free_energy , (float) - 11.1072);
	
	
	// Check bonds
	BOOST_CHECK_EQUAL(segment->bonds[0].first , 0);
	BOOST_CHECK_EQUAL(segment->bonds[0].second , 5);
	
	BOOST_CHECK_EQUAL(segment->bonds[1].first , 1);
	BOOST_CHECK_EQUAL(segment->bonds[1].second , 4);
	
	BOOST_CHECK_EQUAL(segment->bonds[2].first , 2);
	BOOST_CHECK_EQUAL(segment->bonds[2].second , 3);
	
	unlink(filename.c_str());
}

/**
 * @brief checks whether parsing a custom file succeeds - testing both 3' -> 5' && 5' -> 3'
 *
 * @date 2015-06-04
 *
 * @test
 */
BOOST_AUTO_TEST_CASE(Test4)
{
	std::string filename = "tmp.readsegments_test_test4";
	
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile <<   "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
		   "<root>\n"
		   "	<motifs>\n"
		   "		<motif><id>Kt-7 G2nA SAM riboswitch (H. marismortui)</id>\n"
		   "			<sequence_5prime>GAAGAA</sequence_5prime>\n"
		   "			<bonds          >   :::</bonds>\n"
		   "			<sequence_3prime>   AGG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "	</motifs>\n"
		   "</root>\n";
	myfile.close();
	
	
	SegmentTree segments;
	ReadSegments r = ReadSegments(filename, segments);
	
	BOOST_CHECK_EQUAL(segments.empty() , false);
	BOOST_CHECK_EQUAL(segments.size() , 2);
	
	unlink(filename.c_str());
}

/**
 * @brief checks whether parsing a custom file succeeds
 *
 * @date 2015-06-04
 * 
 * @test
 */
BOOST_AUTO_TEST_CASE(Test5)
{
	std::string filename = "tmp.readsegments_test_test5";
	
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile <<   "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
		   "<root>\n"
		   "	<motifs>\n"
		   "		<motif><id>Kt-7 G2nA SAM riboswitch (H. marismortui)</id>\n"
		   "			<sequence_5prime>GAAGAA</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-7 (D. radiodurans)</id>\n"
		   "			<sequence_5prime>UUUGAC</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGC</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-7 (E. coli)</id>\n"
		   "			<sequence_5prime>GUUAUAA</sequence_5prime><bonds          >:    ::</bonds><sequence_3prime>A    AU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-7 (T. thermophilus)</id>\n"
		   "			<sequence_5prime>GUGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-11 (T. thermophilus)</id>\n"
		   "			<sequence_5prime>GACGA C</sequence_5prime><bonds          >   :: :</bonds><sequence_3prime>   ACUA</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-11.eco</id>\n"
		   "			<sequence_5prime>GACGA U</sequence_5prime><bonds          >   :: :</bonds><sequence_3prime>   AUUA</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-15.hma</id>\n"
		   "			<sequence_5prime>AAUG U</sequence_5prime><bonds          >   : :</bonds><sequence_3prime>   AAG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-23.tth</id>\n"
		   "			<sequence_5prime>GCAGAUA</sequence_5prime><bonds          >   ::::</bonds><sequence_3prime>   AUGA</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-23.eco</id>\n"
		   "			<sequence_5prime>GUAGAG</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AUG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-38.hma</id>\n"
		   "			<sequence_5prime>GUUUGAC</sequence_5prime><bonds          >    :::</bonds><sequence_3prime>    AGC</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-42.hma</id>\n"
		   "			<sequence_5prime>CUAGACA</sequence_5prime><bonds          >   :: :</bonds><sequence_3prime>   AG C</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-42.tth</id>\n"
		   "			<sequence_5prime>GAAGACA</sequence_5prime><bonds          >   :: :</bonds><sequence_3prime>   AG C</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-42.dra</id>\n"
		   "			<sequence_5prime>AUAGACA</sequence_5prime><bonds          >   :: :</bonds><sequence_3prime>   AG C</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-42.eco</id>\n"
		   "			<sequence_5prime>CCAGACA</sequence_5prime><bonds          >   :: :</bonds><sequence_3prime>   AG C</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-46.hma</id>\n"
		   "			<sequence_5prime>AUGGAAG</sequence_5prime><bonds          >   ::::</bonds><sequence_3prime>   AGGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-46.tth</id>\n"
		   "			<sequence_5prime>GAUGAAG</sequence_5prime><bonds          >   ::::</bonds><sequence_3prime>   AGGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-46.dra</id>\n"
		   "			<sequence_5prime>GGUGAAG</sequence_5prime><bonds          >   ::::</bonds><sequence_3prime>   AGGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-46.eco</id>\n"
		   "			<sequence_5prime>UGCGAAG</sequence_5prime><bonds          >   ::::</bonds><sequence_3prime>   AGGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-58.hma</id>\n"
		   "			<sequence_5prime>A GGAAG</sequence_5prime><bonds          >   ::::</bonds><sequence_3prime> A AGGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-U4a.hsa</id>\n"
		   "			<sequence_5prime>AAUGA</sequence_5prime><bonds          >   ::</bonds><sequence_3prime>   AG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-U4b.hsa</id>\n"
		   "			<sequence_5prime>AGUGA</sequence_5prime><bonds          >   ::</bonds><sequence_3prime>   AG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.AAG</id>\n"
		   "			<sequence_5prime>AAGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.ACA</id>\n"
		   "			<sequence_5prime>ACAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.ACG</id>\n"
		   "			<sequence_5prime>ACGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.AUA</id>\n"
		   "			<sequence_5prime>AUAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.AUG</id>\n"
		   "			<sequence_5prime>AUGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.AGA</id>\n"
		   "			<sequence_5prime>AGAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.AGG</id>\n"
		   "			<sequence_5prime>AGGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CAA</id>\n"
		   "			<sequence_5prime>CAAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CAG</id>\n"
		   "			<sequence_5prime>CAGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CCA</id>\n"
		   "			<sequence_5prime>CCAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CCG</id>\n"
		   "			<sequence_5prime>CCGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CUA</id>\n"
		   "			<sequence_5prime>CUAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CUG</id>\n"
		   "			<sequence_5prime>CUGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CGU</id><!-- Common -->\n"
		   "			<sequence_5prime>CGUGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UGU</id><!-- Common -->\n"
		   "			<sequence_5prime>UGUGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CGA</id>\n"
		   "			<sequence_5prime>CGAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.CGG</id>\n"
		   "			<sequence_5prime>CGGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UAA</id>\n"
		   "			<sequence_5prime>UAAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UAG</id>\n"
		   "			<sequence_5prime>UAGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UCA</id>\n"
		   "			<sequence_5prime>UCAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UCG</id>\n"
		   "			<sequence_5prime>UCGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UUA</id>\n"
		   "			<sequence_5prime>UUAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UUG</id>\n"
		   "			<sequence_5prime>UUGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UGA</id>\n"
		   "			<sequence_5prime>UGAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.UGG</id>\n"
		   "			<sequence_5prime>UGGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GAA</id>\n"
		   "			<sequence_5prime>GAAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GAG</id>\n"
		   "			<sequence_5prime>GAGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GCA</id>\n"
		   "			<sequence_5prime>GCAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GCG</id>\n"
		   "			<sequence_5prime>GCGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GUA</id>\n"
		   "			<sequence_5prime>GUAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GUG</id>\n"
		   "			<sequence_5prime>GUGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GGA</id>\n"
		   "			<sequence_5prime>GGAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-CD-box.GGG</id>\n"
		   "			<sequence_5prime>GGGGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-L30e.sce</id>\n"
		   "			<sequence_5prime>GCAGAGAU</sequence_5prime><bonds          >+|   ::+</bonds><sequence_3prime>UG   AGG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-SAM-ribo.tte</id>\n"
		   "			<sequence_5prime>GAUGAA</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-T-box.bsu</id>\n"
		   "			<sequence_5prime>AAAGAU</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-c-di-GMP-II.cac</id>\n"
		   "			<sequence_5prime>AAUGAUG</sequence_5prime><bonds          >   :::+</bonds><sequence_3prime>   AGUU</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-G2nA-SAM-riboswitch (T. tengcongensi)</id>\n"
		   "			<sequence_5prime>GACGAA</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AAG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Kt-G2nA-SAMribo.bsu</id>\n"
		   "			<sequence_5prime>GACGAA</sequence_5prime><bonds          >   :::</bonds><sequence_3prime>   AGA</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "		<motif><id>Loop-E-Segment.bac</id>\n"
		   "			<sequence_5prime>GAUGGUA</sequence_5prime><bonds          >:::::::</bonds><sequence_3prime>AUGAGAG</sequence_3prime>\n"
		   "			<energy>-11.1072</energy>\n"
		   "			<directions><five_prime>true</five_prime><three_prime>true</three_prime></directions>\n"
		   "		</motif>\n"
		   "	</motifs>\n"
		   "</root>\n";
	myfile.close();
	
	SegmentTree segment_tree;
	ReadSegments r = ReadSegments(filename, segment_tree);
	
	int n = 61;
	int i = 0;
	
	BOOST_CHECK_EQUAL(segment_tree.empty() , false);
	BOOST_CHECK_EQUAL(segment_tree.size() , n * 2);// forward and reverse
	
	std::pair<Sequence, Sequence> segments[n];
	segments[i++] = {Sequence("GAAGAA"),   Sequence("GGA")};
	segments[i++] = {Sequence("UUUGAC"),   Sequence("CGA")};
	segments[i++] = {Sequence("GUUAUAA"),  Sequence("UAA")};
	segments[i++] = {Sequence("GUGGAU"),   Sequence("GGA")};
	segments[i++] = {Sequence("GACGAC"),   Sequence("AUCA")};
	segments[i++] = {Sequence("GACGAU"),   Sequence("AUUA")};
	segments[i++] = {Sequence("AAUGU"),    Sequence("GAA")};
	segments[i++] = {Sequence("GCAGAUA"),  Sequence("AGUA")};
	segments[i++] = {Sequence("GUAGAG"),   Sequence("GUA")};
	segments[i++] = {Sequence("GUUUGAC"),  Sequence("CGA")};
	segments[i++] = {Sequence("CUAGACA"),  Sequence("CGA")};
	segments[i++] = {Sequence("GAAGACA"),  Sequence("CGA")};
	segments[i++] = {Sequence("AUAGACA"),  Sequence("CGA")};
	segments[i++] = {Sequence("CCAGACA"),  Sequence("CGA")};
	segments[i++] = {Sequence("AUGGAAG"),  Sequence("UGGA")};
	segments[i++] = {Sequence("GAUGAAG"),  Sequence("UGGA")};
	segments[i++] = {Sequence("GGUGAAG"),  Sequence("UGGA")};
	segments[i++] = {Sequence("UGCGAAG"),  Sequence("UGGA")};
	segments[i++] = {Sequence("AGGAAG"),   Sequence("UGGAA")};
	segments[i++] = {Sequence("AAUGA"),    Sequence("GA")};
	segments[i++] = {Sequence("AGUGA"),    Sequence("GA")};
	segments[i++] = {Sequence("AAGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("ACAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("ACGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("AUAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("AUGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("AGAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("AGGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CAAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CAGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CCAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CCGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CUAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CUGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CGUGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UGUGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CGAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("CGGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UAAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UAGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UCAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UCGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UUAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UUGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UGAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("UGGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GAAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GAGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GCAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GCGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GUAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GUGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GGAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GGGGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("GCAGAGAU"), Sequence("GGAGU")};
	segments[i++] = {Sequence("GAUGAA"),   Sequence("GGA")};
	segments[i++] = {Sequence("AAAGAU"),   Sequence("UGA")};
	segments[i++] = {Sequence("AAUGAUG"),  Sequence("UUGA")};
	segments[i++] = {Sequence("GACGAA"),   Sequence("GAA")};
	segments[i++] = {Sequence("GACGAA"),   Sequence("AGA")};
	segments[i++] = {Sequence("GAUGGUA"),  Sequence("GAGAGUA")};
	
	Segment *segment;
	for(i = 0; i < n; i++)
	{
		Position position_5p_s = segments[i].first.data.begin();
		Position position_5p_e = segments[i].first.data.end() - 1;
		Position position_3p_s = segments[i].second.data.begin();
		Position position_3p_e = segments[i].second.data.end() - 1;
		SubSequence subsequence_5p = SubSequence(position_5p_s, position_5p_e);
		SubSequence subsequence_3p = SubSequence(position_3p_s, position_3p_e);
		
		segment = segment_tree.search(subsequence_5p, subsequence_3p);
		
		BOOST_REQUIRE_MESSAGE(segment != NULL , "failed with segment: " << segments[i].first.str().c_str() << "," << segments[i].second.str().c_str());
		BOOST_CHECK_MESSAGE(segment->gibbs_free_energy == (float) - 11.1072 , "failed with segment: " << segments[i].first.str().c_str() << "," << segments[i].second.str().c_str());
	}
	
	unlink(filename.c_str());
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(Testing_Examples)

/**
 * @brief Tests whether there are indeed 28 example sequences within the xml file
 * 
 * @test
 * 
 * @date 2015-06-04
 */
BOOST_AUTO_TEST_CASE(Test_E1)
{
	std::string filename = "share/segmentation-fold/segments.xml";
	SegmentTree segments;
	std::vector<rna_example> rna_examples;
	ReadSegments r = ReadSegments(filename, segments, rna_examples);
	
	BOOST_CHECK_EQUAL(rna_examples.size() , 28);
}

BOOST_AUTO_TEST_SUITE_END()
