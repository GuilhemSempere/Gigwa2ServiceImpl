/*******************************************************************************
 * GIGWA - Service implementation
 * Copyright (C) 2016 - 2019, <CIRAD> <IRD>
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

package fr.cirad.web.filter;

import java.io.IOException;
import java.util.stream.Collectors;

import javax.servlet.FilterChain;
import javax.servlet.FilterConfig;
import javax.servlet.ServletException;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.HttpServletRequest;

import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.service.GigwaGa4ghServiceImpl;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class AssemblyFilter.
 */
public class GigwaAssemblyFilter implements javax.servlet.Filter {

	/** The Constant LOG. */
	private static final Logger LOG = Logger.getLogger(GigwaAssemblyFilter.class);

	/* (non-Javadoc)
	 * @see javax.servlet.Filter#destroy()
	 */
	@Override
	public void destroy()
	{
	}

	/* (non-Javadoc)
	 * @see javax.servlet.Filter#doFilter(javax.servlet.ServletRequest, javax.servlet.ServletResponse, javax.servlet.FilterChain)
	 */
	@Override
	public void doFilter(ServletRequest request, ServletResponse response, FilterChain fc) throws IOException, ServletException
	{
		if (request instanceof HttpServletRequest)
		{
			HttpServletRequest req = (HttpServletRequest) request;
//			Integer assemblyId = null;
			try {
//				String module = req.getParameter("module");
		        String assembly = req.getHeader("assembly");
//		        System.out.println(req.getRequestURI());
//		        System.err.println(request.getReader().lines().collect(Collectors.joining(System.lineSeparator())));
//		        req.getParameterNames().asIterator().forEachRemaining(p -> System.err.println(p));
//		        if (module != null) {
//		        	if (assembly != null) {
//		        		Assembly.setThreadAssembly(MongoTemplateManager.get(module).findOne(new Query(Criteria.where(Assembly.FIELDNAME_NAME).is(assembly)), Assembly.class).getId());
		        		Assembly.setThreadAssembly(assembly == null ? null : Integer.parseInt(assembly));
//		        		LOG.debug(req.getRequestURI() + ": Tied assembly " + Assembly.getThreadAssembly());
//		        	}
//		        	else
//		        		LOG.debug(req.getRequestURI() + ": No assembly identified to set for thread");
//		        }
//		        else
//		        	LOG.debug(req.getRequestURI() + ": No assembly set for thread because no module name");
//				assembly = Integer.parseInt(request.getParameter("assembly"));
			}
			catch (NumberFormatException ignored) {}
//			System.err.println("gigwa " + Assembly.getThreadBoundAssembly());
		}

		fc.doFilter(request, response);
	}

	/* (non-Javadoc)
	 * @see javax.servlet.Filter#init(javax.servlet.FilterConfig)
	 */
	@Override
	public void init(FilterConfig arg0) throws ServletException
	{
	}

}
