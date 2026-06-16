/*******************************************************************************
 * GIGWA - Genotype Investigator for Genome Wide Analyses
 * Copyright (C) 2016 - 2025, <CIRAD> <IRD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 *******************************************************************************/
package fr.cirad.tools.mgdb;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.Collectors;

import fr.cirad.mgdb.model.mongo.maintypes.*;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math.util.MathUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.bson.types.ObjectId;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

//import com.fasterxml.jackson.databind.ObjectMapper;
import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;

import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongo.subtypes.VariantRunDataId;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.model.MgdbSearchVariantsRequest;
import fr.cirad.tools.Helper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class GenotypingDataQueryBuilder.
 */
public class GenotypingDataQueryBuilder implements Iterator<List<BasicDBObject>>
{

    /** The Constant LOG. */
    protected static final Logger LOG = Logger.getLogger(GenotypingDataQueryBuilder.class);

    /** The mongo template. */
    private MongoTemplate mongoTemplate;

    /** The genotyping projects. */
    private Map<Integer, GenotypingProject> genotypingProjects;

    /** Whether or not project has effect annotations. */
    private boolean projectHasEffectAnnotations;

    /** The operator. */
    private List<String> operator;

    private MgdbSearchVariantsRequest req;

    private HashSet<Integer> mutualDiscrim = new HashSet<>();;  // end of the query will be skipped for half of the groups mutually discriminating

    /** The individual or sample to callset list maps (one per group) */
    private List<TreeMap<String /*individual*/, ArrayList<Callset>>> individualOrSampleToCallSetListMaps = new ArrayList<>();

    private int nNextCallCount = 0;

    private int maxAlleleCount = 0;

    private List<Integer> filteredGroups;

    /** Used for shuffling queried chunks */
    private List<Integer> intervalIndexList;

    private List<BasicDBList> intervalQueries = null;

    private BasicDBList variantQueryDBList;

    boolean m_fFilteringOnSequence = false;

    boolean m_fFilteringOnStartSite = false;

    private Document groupFields;

    private boolean fAccountForMultipleRuns = false;

    private boolean fForCounting = false;

    private boolean fGotMultiSampleIndividuals = false;

    private List<String> runsToRestrictQueryTo = null;

    private boolean fExcludeVariantsWithOnlyMissingData = false;

    private Document projectionFields = new Document();

	private Integer m_assemblyId;

    /** The Constant genotypePatternToDescriptionMap. */
    static final private HashMap<String, String> genotypePatternToDescriptionMap = new LinkedHashMap<String, String>();

    /** The Constant genotypePatternToQueryMap. */
    static final private HashMap<String, String> genotypePatternToQueryMap = new HashMap<String, String>();

    static public final String MAIN_RESULT_PROJECTION_FIELD = "r";
    static public final String SECONDARY_RESULT_PROJECTION_FIELD = "r2";

    static
    {
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL, "This will allow all variants whithout applying any filters");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_SAME, "This will allow variants where not all selected individuals have the same genotype");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_MOSTLY_SAME, "This will allow variants where all or most selected individuals have the same genotype");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_DIFFERENT, "This will allow variants where none of the selected individuals have the same genotype");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT, "This will allow variants where some of the selected individuals have the same genotypes");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF, "This will allow variants where selected individuals are all homozygous with the reference allele");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF, "This will allow variants where where at least one selected individual is homozygous with the reference allele");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR, "This will allow variants where selected individuals are all homozygous with an alternate allele");
        genotypePatternToDescriptionMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR, "This will allow variants where at least one selected individual is homozygous with an alternate allele");
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL, null);
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_MOSTLY_SAME, "$eq");
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_SAME, "$eq" + MgdbSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX);
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_DIFFERENT, "$ne");
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_NOT_ALL_DIFFERENT, "$ne" + MgdbSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX);
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_REF, "^0(/0)*$"/*|^$"*/ + MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX);
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_REF, "^0(/0)*$"/*|^$"*/ + MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX);
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ALL_HOMOZYGOUS_VAR, "^([1-9][0-9]*)(/\\1)*$"/*|^$"*/ + MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX);
        genotypePatternToQueryMap.put(MgdbSearchVariantsRequest.GENOTYPE_CODE_LABEL_ATL_ONE_HOMOZYGOUS_VAR, "^([1-9][0-9]*)(/\\1)*$"/*|^$"*/ + MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX);
    }

    public GenotypingDataQueryBuilder(MgdbSearchVariantsRequest gsvr, boolean workWithSamples, BasicDBList variantQueryDBList, boolean fForCounting) throws Exception
    {
        this.req = gsvr;
        mutualDiscrim = new HashSet<>();
        for (int i=0; i<req.getDiscriminate().size(); i++) {
            Integer discrimIndex = req.getDiscriminate().get(i);
            if (discrimIndex != null && i > discrimIndex - 1 && Integer.valueOf(i + 1).equals(req.getDiscriminate().get(discrimIndex - 1)))
                mutualDiscrim.add(i);
        }

        this.fForCounting = fForCounting;
        this.variantQueryDBList = variantQueryDBList;
        Helper.convertIdFiltersToRunFormat(Arrays.asList(this.variantQueryDBList));

    	String info[] = Helper.extractModuleAndProjectIDsFromVariantSetIds(req.getVariantSetId());
        this.m_assemblyId = Assembly.safelyGetThreadBoundAssembly(info[0]);
        for (Object orObject : variantQueryDBList)
            try {
            	List<BasicDBObject> variantFilterAndList = (List<BasicDBObject>) ((BasicDBObject)((Collection)((BasicDBObject) orObject).get("$or")).iterator().next()).get("$and");
            	if (variantFilterAndList != null)		// there can be other $or operators like {"$or": [{"ka": {"$size": 2}}]} for example
	                for (Object variantFilter : variantFilterAndList)
		            	if (((BasicDBObject) variantFilter).containsKey(Assembly.getVariantRefPosPath(m_assemblyId) + "." + ReferencePosition.FIELDNAME_SEQUENCE))
		                    m_fFilteringOnSequence = true;
		                else if (((BasicDBObject) variantFilter).containsKey(Assembly.getVariantRefPosPath(m_assemblyId) + "." + ReferencePosition.FIELDNAME_START_SITE))
		                	m_fFilteringOnStartSite = true;
	        }
	        catch (Exception ignored) {
	        	// probably there isn't such a filter
//	        	ignored.printStackTrace();
	        }

    	List<Integer> projIDs = Arrays.stream(info[1].split(",")).map(pi -> Integer.parseInt(pi)).toList();
        this.mongoTemplate = MongoTemplateManager.get(info[0]);
        this.genotypingProjects = mongoTemplate.find(new Query(Criteria.where("_id").in(projIDs)), GenotypingProject.class).stream().collect(Collectors.toMap(GenotypingProject::getId, Function.identity(), (u, v) -> u, LinkedHashMap::new));

        Query q = new Query();
        q.addCriteria(Criteria.where("_id").in(projIDs));
        q.addCriteria(Criteria.where(GenotypingProject.FIELDNAME_EFFECT_ANNOTATIONS + ".0").exists(true));
        this.projectHasEffectAnnotations = mongoTemplate.findOne(q, GenotypingProject.class) != null;

        int numberGroups = req.getNumberGroups();
//    	if (numberGroups > 0 && genotypingProjects.values().stream().mapToInt(project -> project.getPloidyLevel()).distinct().count() > 1)
//	    	for (int g=0; g<numberGroups; g++)
//	    		if (req.getMaxMaf(g) != null && req.getMaxMaf(g) < 50F || req.getMinMaf(g) != null && req.getMinMaf(g) > 0.0F)
//	    			throw new Exception("MAF filter is only supported across multiple projects if all of them have the same ploidy!");
        this.operator = new ArrayList<>(Collections.nCopies(numberGroups, null));
        this.individualOrSampleToCallSetListMaps = new ArrayList<>(Collections.nCopies(numberGroups, new TreeMap<>()));
        filteredGroups = VariantQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(req, false);
        List<List<String>> callsetIds = req.getAllCallSetIds();
        for (int nGroupIndex = 0; nGroupIndex < callsetIds.size(); nGroupIndex++) {
            this.operator.set(nGroupIndex, genotypePatternToQueryMap.get(req.getGtPattern(nGroupIndex)));
            boolean noMaterialSelected = callsetIds.isEmpty() || callsetIds.get(nGroupIndex).isEmpty();	// will be treated as "all material selected"
            if (workWithSamples) {
                Collection<String> groupSamples = noMaterialSelected ? MgdbDao.getProjectSamples(info[0], projIDs) : callsetIds.get(nGroupIndex).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet());
                this.individualOrSampleToCallSetListMaps.set(nGroupIndex, MgdbDao.getCallsetsBySampleForProjects(info[0], projIDs, groupSamples));
            }
            else {
                Collection<String> groupIndividuals = noMaterialSelected ? MgdbDao.getProjectIndividuals(info[0], projIDs) : callsetIds.get(nGroupIndex).stream().map(csi -> csi.substring(1 + csi.lastIndexOf(Helper.ID_SEPARATOR))).collect(Collectors.toSet());
                this.individualOrSampleToCallSetListMaps.set(nGroupIndex, MgdbDao.getCallsetsByIndividualForProjects(info[0], projIDs, groupIndividuals));
            }
        }

        int nTotalChunkCount = (int) Helper.estimDocCount(mongoTemplate, MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS) + 1;
        if (nTotalChunkCount == 1) {
            MgdbDao.prepareDatabaseForSearches(info[0]);    // list does not exist: create it
            nTotalChunkCount = (int) Helper.estimDocCount(mongoTemplate, MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS) + 1;
        }

        if (m_fFilteringOnSequence && m_fFilteringOnStartSite) {	// when filtering on a precise zone of the genome it's slightly quicker to chunk that given zone rather than use tagged variants
	        Long[] minMaxFound = new Long[2];   // will be filled in by the method below
	        List<String> variantTypes = Arrays.asList(gsvr.getSelectedVariantTypes().split(";"));
	        List<String> sequences = Arrays.asList(gsvr.getReferenceName().split(";"));
	        if (Helper.findDefaultRangeMinMax(info[0], projIDs, null, variantTypes, sequences, gsvr.getStart(), gsvr.getEnd(), minMaxFound))
	        	this.intervalQueries = Helper.getIntervalQueries(nTotalChunkCount, sequences, variantTypes, minMaxFound[0], minMaxFound[1], variantQueryDBList).stream().map(dbo -> { BasicDBList l = new BasicDBList(); l.add(dbo); return l; } ).collect(Collectors.toList());
        }

        if (this.intervalQueries == null) {
        	this.intervalQueries = new ArrayList<>();
        	List<Map> taggedVariantList = mongoTemplate.findAll(Map.class, MgdbDao.COLLECTION_NAME_TAGGED_VARIANT_IDS);
        	for (int i=0; i<nTotalChunkCount; i++) {
        		BasicDBList initialMatchList = new BasicDBList();
        		BasicDBObject idRangeFilter = new BasicDBObject();
        		if (i > 0)
        		    idRangeFilter.append("$gt", taggedVariantList.get(i - 1).get("_id"));
        		if (i < taggedVariantList.size())
        		    idRangeFilter.append("$lte", taggedVariantList.get(i).get("_id"));
        		if (!idRangeFilter.isEmpty())
        			initialMatchList.add(new BasicDBObject("_id", idRangeFilter));
        		if (variantQueryDBList.size() > 0)
        			initialMatchList.addAll(variantQueryDBList);

        		this.intervalQueries.add(initialMatchList);
        	}
        }

        intervalIndexList = new ArrayList<>();
        for (int i=0; i<intervalQueries.size(); i++)
            intervalIndexList.add(i);

        fAccountForMultipleRuns = genotypingProjects.size() > 1 || genotypingProjects.values().iterator().next().getRuns().size() > 1;

        for (TreeMap<String, ArrayList<Callset>> groupIndSamples : individualOrSampleToCallSetListMaps)
            if (groupIndSamples != null && groupIndSamples.values().stream().filter(spList -> spList.size() > 1).findFirst().isPresent()) {
                fGotMultiSampleIndividuals = true;
                break;
            }

        if (!fForCounting || fAccountForMultipleRuns) {
            if (!individualOrSampleToCallSetListMaps.isEmpty()) {
                List<String> involvedSamples = new ArrayList<>();
                for (int filteredGroup : filteredGroups) {
                	HashSet<String> groupSamples = new HashSet<>();
                	for (ArrayList<Callset> callSets : individualOrSampleToCallSetListMaps.get(filteredGroup).values())
                		for (Callset cs : callSets)
                			groupSamples.add(cs.getSampleId());
                	involvedSamples.addAll(groupSamples);
                }

                if (genotypingProjects.size() == 1) {	// this optimization is too complex to handle if several projects are involved
                	int projId = genotypingProjects.keySet().iterator().next();
                    List<String> involvedProjectRuns = Helper.getRunsByProjectFromSampleIDs(info[0], involvedSamples).get(projId);
	                if (involvedProjectRuns.size() < genotypingProjects.get(projId).getRuns().size()) { // not all project runs are involved: adding a filter on the run field will make queries faster
                		runsToRestrictQueryTo = involvedProjectRuns;
	                    fExcludeVariantsWithOnlyMissingData = true;  // some variants may have no data for the selected samples, we don't want to include them
	                }
                }
            }

            String refPosField = m_assemblyId != null ? AbstractVariantData.FIELDNAME_POSITIONS : AbstractVariantData.FIELDNAME_REFERENCE_POSITION;
            if (fAccountForMultipleRuns) {
                groupFields = new Document();
                Document firstContents = new Document();
                firstContents.append(ReferencePosition.FIELDNAME_SEQUENCE, "$" + refPosField + (m_assemblyId != null ? "." + m_assemblyId : "") + "." + ReferencePosition.FIELDNAME_SEQUENCE)
                            .append(ReferencePosition.FIELDNAME_START_SITE, "$" + refPosField + (m_assemblyId != null ? "." + m_assemblyId : "") + "." + ReferencePosition.FIELDNAME_START_SITE)
                            .append(ReferencePosition.FIELDNAME_END_SITE, "$" + refPosField + (m_assemblyId != null ? "." + m_assemblyId : "") + "." + ReferencePosition.FIELDNAME_END_SITE);
                groupFields.put(refPosField, new Document("$first", m_assemblyId != null ? new Document(m_assemblyId.toString(), firstContents) : firstContents));
                groupFields.put(VariantData.FIELDNAME_TYPE, new Document("$first", "$" + VariantData.FIELDNAME_TYPE));
                groupFields.put(VariantData.FIELDNAME_KNOWN_ALLELES, new Document("$first", "$" + VariantData.FIELDNAME_KNOWN_ALLELES));
            }
        }
    }

    /**
     * Gets the number of queries.
     *
     * @return the number of queries
     */
    public int getNumberOfQueries()
    {
        return intervalQueries.size();
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    @Override
    public void remove()
    {
        throw new UnsupportedOperationException("Removal not supported");
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public boolean hasNext()
    {
        return intervalQueries.size() > nNextCallCount;
    }

    public List<Integer> shuffleChunkOrder()
    {
        Collections.shuffle(intervalIndexList);
        return new ArrayList<Integer> (intervalIndexList);
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    /*
     * Changes from old model to new model:
     *
     * 1. No $group stage — one document per variant, all runs in sp[p][r][i]
     * 2. runsToRestrictQueryTo removed from initial match (no longer applicable)
     * 3. Genotype extraction: sp[p][r][i] array lookups instead of sampleGenotypes.callsetId.gt field paths
     * 4. Annotation check: ai[p][r][i] lookups instead of sampleGenotypes.callsetId.additionalInfo.annotation
     *    — same order as old model: annotation check applied per callset BEFORE major genotype resolution
     *    — failing callsets return null, excluded from mode-finding
     * 5. MAF: integer encoding — ref = floor(abs(gt)/100) - abs(gt)%100, tot = floor(abs(gt)/100)
     * 6. Het: gt < 0 replaces regex ^([0-9]+)(\/\1)*$
     * 7. Zygosity filters: integer conditions replace regex on d{g} elements, only for bi-allelic variants
     *    — hom-ref: gt > 0 AND abs(gt) mod 100 == 0
     *    — hom-var: gt > 0 AND abs(gt) mod 100 > 0
     * 8. fAccountForMultipleRuns no longer drives individualSampleGenotypeList or $group
     *    — replaced by fGotMultiSampleIndividuals alone for TOUCHY block
     * 9. fExcludeVariantsWithOnlyMissingData kept — still relevant when subset of runs involved
     */

// ── CONSTRUCTOR CHANGES ───────────────────────────────────────────────────────
// In the constructor, remove the groupFields setup block:
//   if (fAccountForMultipleRuns) { groupFields = new Document(); ... }
// Keep fAccountForMultipleRuns flag itself — still used to determine if individuals
// may have multiple callsets (for TOUCHY block). Also keep fExcludeVariantsWithOnlyMissingData.
// Remove runsToRestrictQueryTo assignment (no longer used in next()).
// ─────────────────────────────────────────────────────────────────────────────

    @Override
    public List<BasicDBObject> next()
    {
        List<BasicDBObject> pipeline = new ArrayList<BasicDBObject>();
        BasicDBList annotationMatchList = new BasicDBList(), finalMatchList = new BasicDBList();
        BasicDBList initialMatchList = intervalQueries.set(intervalIndexList.get(nNextCallCount++).intValue(), null);   // set this one to null to save memory

        /* Step to match variants according to effect annotations */
        if (projectHasEffectAnnotations && (req.getGeneName().length() > 0 || req.getVariantEffect().length() > 0)) {
            if (req.getGeneName().length() > 0) {
                BasicDBObject geneNameDBO;
                if ("-".equals(req.getGeneName()))
                    geneNameDBO = new BasicDBObject("$in", new String[] {"", null});
                else if ("+".equals(req.getGeneName()))
                    geneNameDBO = new BasicDBObject("$regex", "^(?!\\s*$).+");
                else
                    geneNameDBO = new BasicDBObject("$in", Helper.split(req.getGeneName(), ","));
                annotationMatchList.add(new BasicDBObject("vai." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE, geneNameDBO));
            }
            if (req.getVariantEffect().length() > 0)
                annotationMatchList.add(new BasicDBObject("vai." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, new BasicDBObject("$in", Helper.split(req.getVariantEffect(), ","))));
        }

        // Note: runsToRestrictQueryTo removed — new model has one document per variant
        pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", initialMatchList)));

        int maxGroupIndex = filteredGroups.stream().max(Integer::compare).get() + 1;
        boolean[] fZygosityRegex = new boolean[maxGroupIndex];
        boolean[] fNegateMatch = new boolean[maxGroupIndex];
        boolean[] fOr = new boolean[maxGroupIndex];
        boolean[] fMafApplied = new boolean[maxGroupIndex];
        boolean[] fMissingDataApplied = new boolean[maxGroupIndex];
        boolean[] fHezRatioApplied = new boolean[maxGroupIndex];
        boolean[] fCompareBetweenGenotypes = new boolean[maxGroupIndex];
        String[] cleanOperator = new String[maxGroupIndex];

        if (req.getNumberGroups() > 0)
            for (int g : filteredGroups) {
                cleanOperator[g] = operator.get(g);
                if (cleanOperator[g] != null) {
                    if (cleanOperator[g].endsWith(MgdbSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX)) {
                        fNegateMatch[g] = true;
                        cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - MgdbSearchVariantsRequest.AGGREGATION_QUERY_NEGATION_SUFFIX.length());
                    }
                    else if (cleanOperator[g].endsWith(MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX)) {
                        fZygosityRegex[g] = true;
                        cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_ALL_IND_SUFFIX.length());
                    }
                    else if (cleanOperator[g].endsWith(MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX)) {
                        fZygosityRegex[g] = true;
                        fOr[g] = true;
                        cleanOperator[g] = cleanOperator[g].substring(0, cleanOperator[g].length() - MgdbSearchVariantsRequest.AGGREGATION_QUERY_REGEX_APPLY_TO_AT_LEAST_ONE_IND_SUFFIX.length());
                    }
                }

                int nMaxNumberOfAllelesForOneVariant = maxAlleleCount > 0 ? maxAlleleCount : genotypingProjects.values().stream().mapToInt(project -> project.getAlleleCounts().last()).max().getAsInt();
                int nMaxPloidy = genotypingProjects.values().stream().mapToInt(project -> project.getPloidyLevel()).max().getAsInt();
                int nNumberOfPossibleGenotypes;
                try {
                    // nNumberOfPossibleGenotypes = (int) (nMaxNumberOfAllelesForOneVariant + MathUtils.factorial(nMaxNumberOfAllelesForOneVariant) / (MathUtils.factorial(nMaxPloidy) * MathUtils.factorial(nMaxNumberOfAllelesForOneVariant - nMaxPloidy)));
                    nNumberOfPossibleGenotypes = Math.toIntExact((int) MathUtils.factorial(nMaxNumberOfAllelesForOneVariant + nMaxPloidy - 1) / (MathUtils.factorial(nMaxPloidy) * MathUtils.factorial(nMaxNumberOfAllelesForOneVariant - 1)));
                }
                catch (ArithmeticException ae) {    // nMaxNumberOfAllelesForOneVariant must be too large for its factorial to be calculated
                    nNumberOfPossibleGenotypes = Integer.MAX_VALUE;
                }

                int nGroupSize = individualOrSampleToCallSetListMaps.get(g).size();
                double maxMissingGenotypeCount = nGroupSize * req.getMaxMissingData(g) / 100;
                if ("$ne".equals(cleanOperator[g]) && !fNegateMatch[g]) {
                    if (nGroupSize - maxMissingGenotypeCount > nNumberOfPossibleGenotypes) {
                        initialMatchList.add(new BasicDBObject("_id", null));
                        if (nNextCallCount == 1)
                            LOG.info("Aborting 'all different' filter (more called individuals than possible genotypes in group " + (g + 1) + ")");
                        return pipeline;
                    }
                }

                fCompareBetweenGenotypes[g] = cleanOperator[g] != null && !fZygosityRegex[g];
                if ("$ne".equals(cleanOperator[g]) && fNegateMatch[g]) {
                    if (nGroupSize - maxMissingGenotypeCount > nNumberOfPossibleGenotypes) {
                        fCompareBetweenGenotypes[g] = false;
                        if (nNextCallCount == 1)
                            LOG.info("Ignoring 'not all different' filter on group 1 (more called individuals than possible genotypes in group " + (g + 1) + ")");
                    }
                }

                fMafApplied[g] = req.getMaxMaf(g) != null && req.getMaxMaf(g) < 50F || req.getMinMaf(g) != null && req.getMinMaf(g) > 0.0F;
                fMissingDataApplied[g] = (req.getMinMissingData(g) != null && req.getMinMissingData(g) > 0) || (req.getMaxMissingData(g) != null && req.getMaxMissingData(g) < 100);
                fHezRatioApplied[g] = (req.getMinHeZ(g) != null && req.getMinHeZ(g) > 0) || (req.getMaxHeZ(g) != null && req.getMaxHeZ(g) < 100);

                if (req.getMinHeZ(g) != null && req.getMinHeZ(g) > 0 && fZygosityRegex[g]) {
                    if (fOr[g]) {
                        double minMissing = !fMissingDataApplied[g] ? 0 : Math.ceil(req.getMinMissingData(g) * nGroupSize / 100);
                        double minNonHomoZ = Math.ceil((req.getMinHeZ(g)) * (nGroupSize * (1 - (!fMissingDataApplied[g] ? 0 : req.getMinMissingData(g)) / 100f)) / 100);
                        if (minNonHomoZ + minMissing >= nGroupSize) {
                            initialMatchList.add(new BasicDBObject("_id", null));
                            if (nNextCallCount == 1)
                                LOG.info("Aborting incompatible heterozygous / homozygous filters (no space for any homozygous)");
                            return pipeline;
                        }
                    }
                    else {
                        initialMatchList.add(new BasicDBObject("_id", null));
                        if (nNextCallCount == 1)
                            LOG.info("Aborting incompatible heterozygous / homozygous filters (no space for any heterozygous)");
                        return pipeline;
                    }
                }
            }

        BasicDBObject addFieldsVars = new BasicDBObject();
        BasicDBObject addFieldsIn = new BasicDBObject();
        BasicDBObject vars = new BasicDBObject();
        BasicDBObject in = new BasicDBObject();
        BasicDBObject subIn = new BasicDBObject();

        if (req.getNumberGroups() > 0)
            for (int g : filteredGroups) {
                boolean fMostSameSelected = "$eq".equals(cleanOperator[g]) && !fNegateMatch[g];
                boolean fNeedGtArray = fExcludeVariantsWithOnlyMissingData || fMissingDataApplied[g] || fMafApplied[g] || fZygosityRegex[g] || fHezRatioApplied[g] || fCompareBetweenGenotypes[g] || req.isDiscriminate(g);
                int nGroupSize = individualOrSampleToCallSetListMaps.get(g).size();

                List<Object> currentIndGroupGtArray = new ArrayList<>();
                for (String ind : individualOrSampleToCallSetListMaps.get(g).keySet()) {
                    BasicDBList individualSampleGenotypeList = new BasicDBList();
                    List<Callset> individualCallsets = individualOrSampleToCallSetListMaps.get(g).get(ind);
                    if (individualCallsets != null)
                        for (int k = 0; k < individualCallsets.size(); k++) {
                            Callset cs = individualCallsets.get(k);

                            // ── Build sp[p][r][i] lookup ──────────────────────────────────────────
                            int projIdx = cs.getProjectId();
                            List<String> projectRuns = genotypingProjects.get(projIdx).getRuns();
                            int runIdx = projectRuns.indexOf(cs.getRun());
                            int sampleIdx = cs.getCallSetIndexInsideArray();

                            Object spLookup = new BasicDBObject("$arrayElemAt", Arrays.asList(
                                    new BasicDBObject("$arrayElemAt", Arrays.asList(
                                            new BasicDBObject("$arrayElemAt", Arrays.asList(
                                                    "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES, projIdx)),
                                            runIdx)),
                                    sampleIdx));

                            // ── Build annotation conditions from sa[p][r][i] ──────────────────────
                            // Same order as old model: annotation check applied per callset
                            // BEFORE major genotype resolution. Failing → null → treated as missing.
                            BasicDBList annotConditions = new BasicDBList();
                            if (req.getAnnotationFieldThresholds(g) != null) {
                                for (String annotation : req.getAnnotationFieldThresholds(g).keySet()) {
                                    Float threshold = req.getAnnotationFieldThresholds(g).get(annotation);
                                    if (threshold == 0)
                                        continue;
                                    // sa[p][r][i] has same structure as sp[p][r][i]
                                    // but stores {annotationKey: value} objects
                                    Object saLookup = new BasicDBObject("$arrayElemAt", Arrays.asList(
                                            new BasicDBObject("$arrayElemAt", Arrays.asList(
                                                    new BasicDBObject("$arrayElemAt", Arrays.asList(
                                                            "$" + VariantData.SECTION_ADDITIONAL_INFO, projIdx)),
                                                    runIdx)),
                                            sampleIdx));
                                    // Access annotation field: $$annotObj.annotationKey
                                    // saLookup re-evaluated per annotation key — wrap in $let to compute once
                                    annotConditions.add(new BasicDBObject("$lt", Arrays.asList(
                                            "$$annotObj_" + k + "." + annotation, threshold)));
                                }
                            }

                            Object genotypeExpr;
                            if (annotConditions.size() > 0) {
                                // Build sa lookup once via $let, then check all annotation conditions
                                Object saLookup = new BasicDBObject("$arrayElemAt", Arrays.asList(
                                        new BasicDBObject("$arrayElemAt", Arrays.asList(
                                                new BasicDBObject("$arrayElemAt", Arrays.asList(
                                                        "$" + VariantData.SECTION_ADDITIONAL_INFO, projIdx)),
                                                runIdx)),
                                        sampleIdx));
                                // Rebuild annotConditions using $$annotObj reference
                                BasicDBList finalAnnotConditions = new BasicDBList();
                                for (String annotation : req.getAnnotationFieldThresholds(g).keySet()) {
                                    Float threshold = req.getAnnotationFieldThresholds(g).get(annotation);
                                    if (threshold == 0)
                                        continue;
                                    finalAnnotConditions.add(new BasicDBObject("$lt", Arrays.asList(
                                            "$$annotObj." + annotation, threshold)));
                                }
                                genotypeExpr = new BasicDBObject("$let",
                                        new BasicDBObject("vars", new BasicDBObject("annotObj", saLookup))
                                                .append("in", new BasicDBObject("$cond", new Object[] {
                                                        new BasicDBObject("$or", finalAnnotConditions),
                                                        null,
                                                        spLookup
                                                })));
                            }
                            else {
                                genotypeExpr = spLookup;
                            }

                            // fGotMultiSampleIndividuals drives list building (fAccountForMultipleRuns
                            // no longer triggers $group in new model — one doc per variant already)
                            if (fGotMultiSampleIndividuals)
                                individualSampleGenotypeList.add(genotypeExpr);

                            if (k > 0)
                                continue;    // remaining code only once per individual

                            if (fNeedGtArray) {
                                if (fGotMultiSampleIndividuals)
                                    currentIndGroupGtArray.add(individualSampleGenotypeList);
                                else
                                    currentIndGroupGtArray.add(genotypeExpr);
                            }
                        }
                }

                // ── MAF — integer encoding ────────────────────────────────────────────────
                // ref = floor(abs(gt)/100) - (abs(gt) mod 100)
                // tot = floor(abs(gt)/100)
                if (fMafApplied[g]) {
                    BasicDBObject reduceIn = new BasicDBObject()
                            .append("ref", new BasicDBObject("$add", Arrays.asList(
                                    "$$value.ref",
                                    new BasicDBObject("$subtract", Arrays.asList(
                                            new BasicDBObject("$floor", new BasicDBObject("$divide", Arrays.asList(
                                                    new BasicDBObject("$abs", "$$this"), 100))),
                                            new BasicDBObject("$mod", Arrays.asList(
                                                    new BasicDBObject("$abs", "$$this"), 100)))))))
                            .append("tot", new BasicDBObject("$add", Arrays.asList(
                                    "$$value.tot",
                                    new BasicDBObject("$floor", new BasicDBObject("$divide", Arrays.asList(
                                            new BasicDBObject("$abs", "$$this"), 100))))));

                    in.put("ac" + g, new BasicDBObject("$reduce", new BasicDBObject()
                            .append("input", "$$gt" + g)
                            .append("initialValue", new BasicDBObject("ref", 0).append("tot", 0))
                            .append("in", reduceIn)));
                }

                // ── Missing count — unchanged ─────────────────────────────────────────────
                if (fExcludeVariantsWithOnlyMissingData || fMissingDataApplied[g] || fHezRatioApplied[g] || (cleanOperator[g] != null && !fZygosityRegex[g]) || req.isDiscriminate(g)) {
                    in.put("m" + g, new BasicDBObject("$subtract", Arrays.asList(nGroupSize, new BasicDBObject("$size", "$$gt" + g))));
                }

                // ── Distinct count — unchanged (integer equality same as string) ──────────
                if (fCompareBetweenGenotypes[g] && !fMostSameSelected) {
                    BasicDBObject filter = new BasicDBObject("input", new BasicDBObject("$setUnion", "$$gt" + g));
                    filter.put("as", "gt");
                    filter.put("cond", new BasicDBObject("$ne", Arrays.asList("$$gt", null)));
                    in.put("dc" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter)));
                }

                // ── Distinct set — unchanged ──────────────────────────────────────────────
                if (fZygosityRegex[g] || fMostSameSelected || req.isDiscriminate(g)) {
                    if (fMostSameSelected || req.isDiscriminate(g))
                        in.put("gt" + g, "$$gt" + g);
                    in.put("d" + g, new BasicDBObject("$setUnion", "$$gt" + g));
                }

                // ── Missing data filter — unchanged ───────────────────────────────────────
                if (fMissingDataApplied[g]) {
                    BasicDBObject missingDataFilter = new BasicDBObject();
                    if (req.getMinMissingData(g) != null && req.getMinMissingData(g) > 0)
                        missingDataFilter.put("$gte", nGroupSize * req.getMinMissingData(g) / 100);
                    if (req.getMaxMissingData(g) != null && req.getMaxMissingData(g) < 100)
                        missingDataFilter.put("$lte", nGroupSize * req.getMaxMissingData(g) / 100);
                    finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".m" + g, missingDataFilter));
                }

                // ── Het ratio — gt < 0 replaces regex ────────────────────────────────────
                if (fHezRatioApplied[g]) {
                    BasicDBObject filter = new BasicDBObject("input", "$$gt" + g);
                    filter.put("as", "gt");
                    // heterozygous = negative integer in new encoding
                    filter.put("cond", new BasicDBObject("$lt", Arrays.asList("$$gt", 0)));
                    in.put("he" + g, new BasicDBObject("$size", new BasicDBObject("$filter", filter)));

                    BasicDBObject hzFilter = new BasicDBObject();
                    if (req.getMinHeZ(g) != null && req.getMinHeZ(g) > 0)
                        hzFilter.put("$gte", req.getMinHeZ(g));
                    if (req.getMaxHeZ(g) != null && req.getMaxHeZ(g) < 100)
                        hzFilter.put("$lte", req.getMaxHeZ(g));
                    finalMatchList.add(new BasicDBObject(SECONDARY_RESULT_PROJECTION_FIELD + ".hef" + g, hzFilter));
                }

                boolean innerLetNeeded = false;
                for (boolean[] booleans : Arrays.asList(new boolean[] {fExcludeVariantsWithOnlyMissingData, fMostSameSelected, req.getDiscriminate().stream().anyMatch(n -> n != null)}, fZygosityRegex, fMafApplied, fCompareBetweenGenotypes, fHezRatioApplied))
                    if (ArrayUtils.contains(booleans, true)) {
                        innerLetNeeded = true;
                        break;
                    }

                if (innerLetNeeded) {
                    // keep previously computed fields
                    if (fExcludeVariantsWithOnlyMissingData || fMissingDataApplied[g] || fCompareBetweenGenotypes[g] || fHezRatioApplied[g] || req.isDiscriminate(g))
                        subIn.put("m" + g, "$$m" + g);
                    if (fZygosityRegex[g] || fMostSameSelected || req.isDiscriminate(g))
                        subIn.put("d" + g, "$$d" + g);
                    if (fCompareBetweenGenotypes[g] && !fMostSameSelected)
                        subIn.put("dc" + g, "$$dc" + g);
                    if (fHezRatioApplied[g])
                        subIn.put("he" + g, "$$he" + g);

                    if (fCompareBetweenGenotypes[g] && !fMostSameSelected) {
                        subIn.put("dm" + g, new BasicDBObject("$add", new Object[] {"$$dc" + g, "$$m" + g}));

                        finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".m" + g, new BasicDBObject("$lt", nGroupSize - 1)));
                        if ("$eq".equals(cleanOperator[g]) && fNegateMatch[g])
                            finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".dc" + g, new BasicDBObject("$ne", 1)));
                        else if ("$ne".equals(cleanOperator[g]))
                            finalMatchList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".dm" + g, new BasicDBObject(fNegateMatch[g] ? "$lt" : "$eq", nGroupSize)));
                        else
                            LOG.error("Invalid operator: " + operator);
                    }

                    // ── MAF frequency — same formula, different inputs ────────────────────
                    if (fMafApplied[g]) {
                        BasicDBObject condFrq = new BasicDBObject("$cond", new Object[] {
                                new BasicDBObject("$eq", new Object[] {"$$ac" + g + ".tot", 0}),
                                null,
                                new BasicDBObject("$multiply", new Object[] {
                                        new BasicDBObject("$divide", new Object[] {"$$ac" + g + ".ref", "$$ac" + g + ".tot"}),
                                        100})});
                        subIn.put("f" + g, condFrq);

                        // MAF filter logic unchanged — same OR pattern for ref/alt symmetry
                        float minMaf = req.getMinMaf(g), maxMaf = req.getMaxMaf(g);
                        BasicDBList orMafMatch = new BasicDBList();
                        BasicDBList andMafMatch = new BasicDBList();
                        if (minMaf > 0 && minMaf != maxMaf)
                            andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject("$gte", minMaf)));
                        andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject(minMaf == maxMaf ? "$eq" : "$lte", maxMaf != 50 ? maxMaf : (100F - minMaf))));
                        orMafMatch.add(andMafMatch.size() == 1 ? andMafMatch.iterator().next() : new BasicDBObject("$and", andMafMatch));
                        if (maxMaf != 50) {
                            andMafMatch = new BasicDBList();
                            if (minMaf > 0 && minMaf != maxMaf)
                                andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject("$lte", (100F - minMaf))));
                            andMafMatch.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".f" + g, new BasicDBObject(minMaf == maxMaf ? "$eq" : "$gte", (100F - maxMaf))));
                            orMafMatch.add(andMafMatch.size() == 1 ? andMafMatch.iterator().next() : new BasicDBObject("$and", andMafMatch));
                        }
                        finalMatchList.add(new BasicDBObject("$or", orMafMatch));
                    }
                }

                if (cleanOperator[g] != null || req.isDiscriminate(g)) {
                    if (nGroupSize >= 1) {
                        if (fZygosityRegex[g]) {
                            // ── Zygosity — integer conditions replace regex ────────────────────
                            // cleanOperator[g] is "^0(/0)*$" (hom-ref) or "^([1-9][0-9]*)(/\\1)*$" (hom-var)
                            boolean isHomRef = "^0(/0)*$".equals(cleanOperator[g]);

                            if (!fOr[g]) {
                                // All individuals must match: d{g} has exactly one element AND it satisfies condition
                                // Replaces: d{g}[0] matches regex AND d{g}[1] doesn't exist
                                List<Object> allMatchConditions = new ArrayList<>();
                                // exactly one distinct genotype
                                allMatchConditions.add(new BasicDBObject("$expr",
                                        new BasicDBObject("$eq", Arrays.asList(
                                                new BasicDBObject("$size", "$" + MAIN_RESULT_PROJECTION_FIELD + ".d" + g), 1))));
                                // that genotype satisfies the zygosity condition
                                allMatchConditions.add(new BasicDBObject("$expr",
                                        buildZygosityCondOnArrayElem("$" + MAIN_RESULT_PROJECTION_FIELD + ".d" + g, 0, isHomRef)));
                                finalMatchList.add(new BasicDBObject("$and", allMatchConditions));
                            }
                            else {
                                // At least one individual matches: any element in d{g} satisfies condition
                                // Replaces: OR of d{g}[j] matches regex for j in 0..nGroupSize-1
                                finalMatchList.add(new BasicDBObject("$expr",
                                        new BasicDBObject("$gt", Arrays.asList(
                                                new BasicDBObject("$size", new BasicDBObject("$filter", new BasicDBObject()
                                                        .append("input", "$" + MAIN_RESULT_PROJECTION_FIELD + ".d" + g)
                                                        .append("as", "gt")
                                                        .append("cond", buildZygosityIntegerCond(isHomRef)))),
                                                0))));
                            }
                        }

                        if (fMostSameSelected || req.isDiscriminate(g)) {
                            // c{g} = per-distinct-gt count — unchanged (integer equality same as string)
                            BasicDBObject filter = new BasicDBObject("input", "$$gt" + g);
                            filter.put("as", "g");
                            filter.put("cond", new BasicDBObject("$eq", Arrays.asList("$$g", "$$d")));
                            subIn.put("c" + g, new BasicDBObject("$map",
                                    new BasicDBObject("input", "$$d" + g)
                                            .append("as", "d")
                                            .append("in", new BasicDBObject("$size", new BasicDBObject("$filter", filter)))));

                            addFieldsVars.put("dgc" + g, new BasicDBObject("$max", "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + g));
                            Object minimumDominantGenotypeCount = new BasicDBObject("$multiply", Arrays.asList(
                                    new BasicDBObject("$subtract", new Object[] {nGroupSize, "$" + MAIN_RESULT_PROJECTION_FIELD + ".m" + g}),
                                    req.getMostSameRatio(g) / 100f));

                            if (fMostSameSelected)
                                addFieldsIn.put("ed" + g, new BasicDBObject("$gte", Arrays.asList("$$dgc" + g, minimumDominantGenotypeCount)));

                            if (!mutualDiscrim.contains(g)) {
                                Integer groupToDiscriminateWith = req.getDiscriminate().get(g);
                                if (groupToDiscriminateWith != null) {
                                    groupToDiscriminateWith--;
                                    // Discriminate logic unchanged — $ne on integers same as on strings
                                    addFieldsIn.put("dd" + g, new BasicDBObject("$and", Arrays.asList(
                                            new BasicDBObject("$eq", Arrays.asList(1, new BasicDBObject("$size", new BasicDBObject("$filter",
                                                    new BasicDBObject("input", "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + groupToDiscriminateWith)
                                                            .append("cond", new BasicDBObject("$eq", Arrays.asList("$$dgc" + groupToDiscriminateWith, "$$this"))))))),
                                            new BasicDBObject("$eq", Arrays.asList(1, new BasicDBObject("$size", new BasicDBObject("$filter",
                                                    new BasicDBObject("input", "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + g)
                                                            .append("cond", new BasicDBObject("$eq", Arrays.asList("$$dgc" + g, "$$this"))))))),
                                            new BasicDBObject("$ne", Arrays.asList(
                                                    new BasicDBObject("$arrayElemAt", Arrays.asList(
                                                            "$" + MAIN_RESULT_PROJECTION_FIELD + ".d" + groupToDiscriminateWith,
                                                            new BasicDBObject("$indexOfArray", Arrays.asList(
                                                                    "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + groupToDiscriminateWith,
                                                                    "$$dgc" + groupToDiscriminateWith)))),
                                                    new BasicDBObject("$arrayElemAt", Arrays.asList(
                                                            "$" + MAIN_RESULT_PROJECTION_FIELD + ".d" + g,
                                                            new BasicDBObject("$indexOfArray", Arrays.asList(
                                                                    "$" + MAIN_RESULT_PROJECTION_FIELD + ".c" + g,
                                                                    "$$dgc" + g))))))
                                    )));
                                    finalMatchList.add(new BasicDBObject(SECONDARY_RESULT_PROJECTION_FIELD + ".dd" + g, true));
                                }
                            }

                            if (fMostSameSelected)
                                finalMatchList.add(new BasicDBObject(SECONDARY_RESULT_PROJECTION_FIELD + ".ed" + g, true));
                        }
                    }
                }

                if (fHezRatioApplied[g]) {
                    BasicDBList condList = new BasicDBList();
                    condList.add(new BasicDBObject("$eq", new Object[] {"$" + MAIN_RESULT_PROJECTION_FIELD + ".m" + g, nGroupSize}));
                    condList.add(0);
                    condList.add(new BasicDBObject("$divide", Arrays.asList(
                            new BasicDBObject("$multiply", new Object[] {"$" + MAIN_RESULT_PROJECTION_FIELD + ".he" + g, 100}),
                            new BasicDBObject("$subtract", new Object[] {nGroupSize, "$" + MAIN_RESULT_PROJECTION_FIELD + ".m" + g}))));
                    addFieldsIn.put("hef" + g, new BasicDBObject("$cond", condList));
                }

                // ── TOUCHY block — unchanged logic, works with integers ───────────────────────
                // Mode-finding uses value equality ($eq) which works the same for integers
                Object groupGenotypesArray;
                if (fNeedGtArray && fGotMultiSampleIndividuals) {

                    // Step 1: for each distinct value in $$ig, count non-null occurrences
                    //         produces [{k: genotypeValue, v: count}, ...]
                    BasicDBObject countFilter = new BasicDBObject("input", "$$ig")
                            .append("as", "el")
                            .append("cond", new BasicDBObject("$and", Arrays.asList(
                                    new BasicDBObject("$ne", Arrays.asList("$$el", null)),
                                    new BasicDBObject("$eq", Arrays.asList("$$el", "$$it")))));

                    BasicDBObject countsArray = new BasicDBObject("$map",
                            new BasicDBObject("input", new BasicDBObject("$setUnion", Arrays.asList("$$ig")))
                                    .append("as", "it")
                                    .append("in", new BasicDBObject("k", "$$it")
                                            .append("v", new BasicDBObject("$size",
                                                    new BasicDBObject("$filter", countFilter)))));

                    // Step 2: keep only entries whose count equals max count (major gt candidates)
                    BasicDBObject maxCountItems = new BasicDBObject("$filter",
                            new BasicDBObject("input", "$$ca")
                                    .append("as", "c")
                                    .append("cond", new BasicDBObject("$eq", Arrays.asList(
                                            "$$c.v", new BasicDBObject("$max", "$$ca.v")))));

                    // Step 3: if tie (>1 winner) → null, else take the winning key
                    BasicDBObject pickMajor = new BasicDBObject("$cond", Arrays.asList(
                            new BasicDBObject("$gt", Arrays.asList(
                                    new BasicDBObject("$size", "$$mci"), 1)),
                            null,
                            new BasicDBObject("$arrayElemAt", Arrays.asList("$$mci.k", 0))));

                    // Step 4: inner $let — bind ca, resolve major gt or null if ig is empty
                    BasicDBObject innerLet = new BasicDBObject("$let",
                            new BasicDBObject("vars", new BasicDBObject("ca", countsArray))
                                    .append("in", new BasicDBObject("$cond", Arrays.asList(
                                            new BasicDBObject("$eq", Arrays.asList(
                                                    new BasicDBObject("$size", "$$ig"), 0)),
                                            null,
                                            new BasicDBObject("$let",
                                                    new BasicDBObject("vars", new BasicDBObject("mci", maxCountItems))
                                                            .append("in", pickMajor))))));

                    // Step 5: map over each individual's observation list to resolve major gt
                    BasicDBObject mapIndividuals = new BasicDBObject("$map",
                            new BasicDBObject("input", "$$igl")
                                    .append("as", "ig")
                                    .append("in", innerLet));

                    // Step 6: outer $let — bind igl = currentIndGroupGtArray, map over individuals
                    groupGenotypesArray = new BasicDBObject("$let",
                            new BasicDBObject("vars", new BasicDBObject("igl", currentIndGroupGtArray))
                                    .append("in", mapIndividuals));
                }
                else {
                    groupGenotypesArray = currentIndGroupGtArray;
                }

                // filter out null genotypes (missing + annotation-failed)
                vars.put("gt" + g, new BasicDBObject("$filter", new BasicDBObject()
                        .append("input", groupGenotypesArray)
                        .append("as", "g")
                        .append("cond", new BasicDBObject("$ne", new Object[] {"$$g", null}))));
            }

        if (subIn.size() > 0) {
            BasicDBObject subVars = in;
            BasicDBObject subLet = new BasicDBObject("vars", subVars);
            subLet.put("in", subIn);
            in = new BasicDBObject("$let", subLet);
        }

        if (in.size() > 0) {
            BasicDBObject let = new BasicDBObject("vars", vars);
            let.put("in", in);
            projectionFields.put(MAIN_RESULT_PROJECTION_FIELD, new BasicDBObject("$let", let));
        }

        // Note: no $group stage in new model — one document per variant already
        if (!fForCounting)
            projectionFields.put("_id", "$_id");

        if (!annotationMatchList.isEmpty())
            pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", annotationMatchList)));

        if (!projectionFields.isEmpty()) {
            projectionFields.put(Assembly.getVariantRefPosPath(m_assemblyId), 1);
            projectionFields.put(AbstractVariantData.FIELDNAME_TYPE, 1);
            projectionFields.put(AbstractVariantData.FIELDNAME_KNOWN_ALLELES, 1);
            pipeline.add(new BasicDBObject("$project", projectionFields));
        }

        if (addFieldsIn.size() > 0) {
            BasicDBObject addFieldsLet = new BasicDBObject("vars", addFieldsVars);
            addFieldsLet.put("in", addFieldsIn);
            pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(SECONDARY_RESULT_PROJECTION_FIELD, new BasicDBObject("$let", addFieldsLet))));
        }

        if (fExcludeVariantsWithOnlyMissingData) {
            ArrayList<BasicDBObject> orList = new ArrayList<BasicDBObject>();
            for (int filteredGroup : filteredGroups) {
                if (req.getMaxMissingData(filteredGroup) == 100)
                    orList.add(new BasicDBObject(MAIN_RESULT_PROJECTION_FIELD + ".m" + filteredGroup,
                            new BasicDBObject("$lt", individualOrSampleToCallSetListMaps.get(filteredGroup).size())));
            }
            if (!orList.isEmpty())
                finalMatchList.add(new BasicDBObject("$or", orList));
        }

        if (finalMatchList.size() > 0)
            pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", finalMatchList)));
        LOG.debug("filtering with: " + pipeline);
        return pipeline;
    }

// ── Helper methods ────────────────────────────────────────────────────────────

    /**
     * Builds an integer zygosity condition to use inside $filter ($$gt variable).
     * Replaces regex matching on string genotypes.
     * Only valid for biallelic variants (incremental multiallelic encoding
     * does not store ploidy/alt count).
     *
     * @param isHomRef true for hom-ref (^0(/0)*$), false for hom-var (^([1-9][0-9]*)(/\1)*$)
     */
    private BasicDBObject buildZygosityIntegerCond(boolean isHomRef) {
        if (isHomRef) {
            // hom-ref: positive (homozygous) AND alt count == 0
            return new BasicDBObject("$and", Arrays.asList(
                    new BasicDBObject("$gt", Arrays.asList("$$gt", 0)),
                    new BasicDBObject("$eq", Arrays.asList(
                            new BasicDBObject("$mod", Arrays.asList(new BasicDBObject("$abs", "$$gt"), 100)), 0))));
        }
        else {
            // hom-var: positive (homozygous) AND alt count > 0
            return new BasicDBObject("$and", Arrays.asList(
                    new BasicDBObject("$gt", Arrays.asList("$$gt", 0)),
                    new BasicDBObject("$gt", Arrays.asList(
                            new BasicDBObject("$mod", Arrays.asList(new BasicDBObject("$abs", "$$gt"), 100)), 0))));
        }
    }

    /**
     * Builds an integer zygosity condition on a specific array element (for $expr in $match).
     * Replaces d{g}[index] regex match.
     *
     * @param arrayPath  MongoDB path to the array (e.g. "$r.d0")
     * @param index      element index
     * @param isHomRef   true for hom-ref, false for hom-var
     */
    private BasicDBObject buildZygosityCondOnArrayElem(String arrayPath, int index, boolean isHomRef) {
        Object elem = new BasicDBObject("$arrayElemAt", Arrays.asList(arrayPath, index));
        if (isHomRef) {
            return new BasicDBObject("$and", Arrays.asList(
                    new BasicDBObject("$gt", Arrays.asList(elem, 0)),
                    new BasicDBObject("$eq", Arrays.asList(
                            new BasicDBObject("$mod", Arrays.asList(new BasicDBObject("$abs", elem), 100)), 0))));
        }
        else {
            return new BasicDBObject("$and", Arrays.asList(
                    new BasicDBObject("$gt", Arrays.asList(elem, 0)),
                    new BasicDBObject("$gt", Arrays.asList(
                            new BasicDBObject("$mod", Arrays.asList(new BasicDBObject("$abs", elem), 100)), 0))));
        }
    }
//    private List<Comparable> buildFullListFromRange(MongoTemplate mongoTemplate, Comparable leftBound, Comparable rightBound) {
//        if (leftBound == null)
//            leftBound = mongoTemplate.findOne(new Query().with(new Sort(Sort.Direction.ASC, "_id")), VariantData.class).getId();
//        if (rightBound == null)
//            rightBound = mongoTemplate.findOne(new Query().with(new Sort(Sort.Direction.DESC, "_id")), VariantData.class).getId();
//
//        ArrayList<Comparable> result = new ArrayList<>();
//        String leftAsString = leftBound.toString(), rightAsString = rightBound.toString();
//        if (ObjectId.isValid(leftAsString) && ObjectId.isValid(rightAsString) && leftAsString.substring(0, 18).equals(rightAsString.substring(0, 18)))
//        {
//            int nCurrentId = Integer.parseInt(leftAsString.substring(18, 24), 16);
//            while (nCurrentId <= Integer.parseInt(rightAsString.substring(18, 24), 16))
//                result.add(new ObjectId(leftAsString.substring(0, 18) + Integer.toHexString(nCurrentId++)));
//        }
//        else
//        {
//            Query q = new Query().with(new Sort(Sort.Direction.ASC, "_id"));
//            q.fields().include("_id");
//            // not finished implementing method (functionality non needed)
//        }
//        return result;
//    }

    public static boolean areObjectIDsConsecutive(ObjectId first, ObjectId second)
    {
//        if (first == null || second == null)
//            return false;

        String firstAsString = first.toHexString(), secondAsString = second.toHexString();
        if (!firstAsString.substring(0, 18).equals(secondAsString.substring(0, 18)))
            return false;

        return 1 + Integer.parseInt(firstAsString.substring(18, 24), 16) == Integer.parseInt(secondAsString.substring(18, 24), 16);
    }

//    public static Document tryAndShrinkIdList(String pathToVariantId, Collection<Comparable> idCollection, int nShrinkThreshold)
//    {
//        if (idCollection.size() >= 300000)
//            try
//            {
//        //        long b4 = System.currentTimeMillis();
//        //        SortedSet<Comparable> idSet = SortedSet.class.isAssignableFrom(idCollection.getClass()) ? (SortedSet<Comparable>) idCollection : new TreeSet<Comparable>(idCollection);
//        //        System.out.println("sorting took " + (System.currentTimeMillis() - b4));
//                BasicDBList orList = new BasicDBList();
//                ArrayList<ObjectId> inIdList = new ArrayList<>(), rangeIdList = new ArrayList<>();
//
//                ObjectId previousId = null;
//                for (Comparable id : idCollection)
//                {
//                    ObjectId currentId = (ObjectId) id;
//                    if (previousId == null || areObjectIDsConsecutive(previousId, currentId))
//                        rangeIdList.add(currentId);
//                    else
//                    {
//                        if (rangeIdList.size() >= nShrinkThreshold)
//                        {    // replace list with a range
//                            BasicDBList chunkMatchAndList = new BasicDBList();
//                            chunkMatchAndList.add(new Document(pathToVariantId, new Document("$gte", rangeIdList.get(0))));
//                            chunkMatchAndList.add(new Document(pathToVariantId, new Document("$lte", rangeIdList.get(rangeIdList.size() - 1))));
//                            orList.add(new Document("$and", chunkMatchAndList));
//                        }
//                        else
//                            inIdList.addAll(rangeIdList);    // range is too small, keep the list
//
//                        rangeIdList.clear();
//                        rangeIdList.add(currentId);
//                    }
//                    previousId = currentId;
//                }
//                inIdList.addAll(rangeIdList);
//
//                if (inIdList.size() > 0 || orList.size() == 0)
//                    orList.add(new Document(pathToVariantId, new Document("$in", inIdList)));
//
//                return orList.size() > 1 ? new Document("$or", orList) : (Document) orList.iterator().next();
//            }
//            catch (ClassCastException cce)
//            {
//                if (!cce.getMessage().contains("ObjectId"))
//                    throw cce;    // otherwise it simply means IDs are of a different type, in which case we can't shrink the collection
//            }
////        else
////        {
////            LOG.debug("Didn't shrink id collection (" + idCollection.size() + " records only)");
////        }
//
//        return new Document(pathToVariantId, new Document("$in", idCollection));    // not shrinked
//    }

    public static HashMap<String, String> getGenotypePatternToQueryMap() {
        return genotypePatternToQueryMap;
    }

    public static HashMap<String, String> getGenotypePatternToDescriptionMap() {
        return genotypePatternToDescriptionMap;
    }

    public void setMaxAlleleCount(int maxAlleleCount) {
        this.maxAlleleCount = maxAlleleCount;
    }
}