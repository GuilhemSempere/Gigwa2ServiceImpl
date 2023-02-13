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
public class AssemblyFilter implements javax.servlet.Filter {

	/** The Constant LOG. */
	private static final Logger LOG = Logger.getLogger(AssemblyFilter.class);

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

			Integer assemblyId = null;
			try {
		        String assembly = ((HttpServletRequest) request).getHeader("assembly");
		        assemblyId = assembly == null ? null : MongoTemplateManager.get(request.getParameter("module")).findOne(new Query(Criteria.where(Assembly.FIELDNAME_NAME).is(assembly)), Assembly.class).getId();
//				assembly = Integer.parseInt(request.getParameter("assembly"));
			}
			catch (NumberFormatException ignored) {}
			GigwaGa4ghServiceImpl.setThreadAssembly(assemblyId);
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
