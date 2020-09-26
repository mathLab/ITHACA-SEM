<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow 3DH1D P=2 Boundary Conditions from files</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_bcsfromfiles.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH1D_bcsfromfiles.xml</file>
        <file description="Session File">ChanFlow_3DH1D_bcsfromfiles_walls.bc</file>
	<file description="Session File">ChanFlow_3DH1D_bcsfromfiles_inflow.bc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">1.40307e-16</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>            
	    <value variable="p" tolerance="1e-6">6.68056e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">1.249e-15</value>
            <value variable="v" tolerance="1e-6">1.89879e-16</value>
            <value variable="w" tolerance="1e-6">9.51393e-18</value>            
	    <value variable="p" tolerance="1e-6">3.68594e-14</value>
        </metric>
    </metrics>
</test>


