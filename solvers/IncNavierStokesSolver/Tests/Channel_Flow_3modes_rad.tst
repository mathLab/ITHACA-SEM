<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow 2D with Radiation outflow </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.18976e-16</value>
            <value variable="v" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="1e-12">3.23754e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.02696e-15</value>
            <value variable="v" tolerance="1e-12">5.72397e-17</value>
	    <value variable="p" tolerance="1e-12">4.21885e-15</value>
        </metric>
    </metrics>
</test>


