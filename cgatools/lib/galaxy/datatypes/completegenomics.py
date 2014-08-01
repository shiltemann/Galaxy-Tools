"""
Complete Genomics datatypes
Birgit Crain - Complete Genomics, Inc
"""

import pkg_resources
pkg_resources.require( "bx-python" )

import logging
from galaxy.datatypes import data
from galaxy import util
from cgi import escape
from galaxy.datatypes import metadata
from galaxy.datatypes import tabular
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.tabular import Tabular
import galaxy_utils.sequence.vcf
from galaxy.datatypes.sniff import *

log = logging.getLogger(__name__)

class CG_Var( Tabular ):
    file_ext = 'cg_var'
    def __init__(self, **kwd):
        """Initialize CG_Var datatype"""
        Tabular.__init__( self, **kwd )
        self.column_names = ['locus', 'ploidy', 'allele', 'chromosome', 'begin', 'end',
                             'varType', 'reference', 'alleleSeq', 'varScoreVAF',
                             'varScoreEAF', 'varQuality', 'hapLink', 'xRef'
                             ]
    def display_peek( self, dataset ):
        """Returns formated html of peek"""
        return Tabular.make_html_table( self, dataset, column_names=self.column_names )

class CG_MasterVar( Tabular ):
    file_ext = 'cg_mastervar'
    def __init__(self, **kwd):
        """Initialize CG_MasterVar datatype"""
        Tabular.__init__( self, **kwd )
        self.column_names = ['locus', 'ploidy', 'chromosome', 'begin', 'end', 'zygosity',
                             'varType', 'reference', 'allele1Seq', 'allele2Seq',
                             'allele1VarScoreVAF', 'allele2VarScoreVAF', 'allele1VarScoreEAF',
                             'allele2VarScoreEAF', 'allele1VarQuality', 'allele2VarQuality',
                             'allele1HapLink', 'allele2HapLink', 'allele1XRef', 'allele2XRef',
                             'evidenceIntervalId', 'allele1ReadCount', 'allele2ReadCount',
                             'referenceAlleleRead', 'totalReadCount', 'allele1Gene',
                             'allele2Gene	pfam', 'miRBaseId', 'repeatMasker', 'segDupOverlap',
                             'relativeCoverageDiploid', 'calledPloidy',
                             'relativeCoverageNondiploid', 'calledLevel'
                             ]
    
    def display_peek( self, dataset ):
        """Returns formated html of peek"""
        return Tabular.make_html_table( self, dataset, column_names=self.column_names )
        
class CG_Gene( Tabular ):
    file_ext = 'cg_gene'
    def __init__(self, **kwd):
        """Initialize CG_Gene datatype"""
        Tabular.__init__( self, **kwd )
        self.column_names = ['index', 'locus', 'allele', 'chromosome', 'begin', 'end',
                             'varType', 'reference', 'call', 'xRef', 'geneId',
                             'mrnaAcc', 'proteinAcc', 'symbol', 'orientation', 'component',
                             'componentIndex', 'hasCodingRegion', 'impact', 'nucleotidePos',
                             'proteinPos', 'annotationRefSequence', 'sampleSequence',
                             'genomeRefSequence', 'pfam'
                             ]
   
    def display_peek( self, dataset ):
        """Returns formated html of peek"""
        return Tabular.make_html_table( self, dataset, column_names=self.column_names )
