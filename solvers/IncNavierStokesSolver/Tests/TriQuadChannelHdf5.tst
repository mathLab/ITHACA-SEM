<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Tri and Quad element channel Flow 2D reading hdf5 (3 proc) </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>TriQuadChannelHdf5.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">TriQuadChannelHdf5.xml</file>
        <file description="Restart File">TriQuadChannelHdf5.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">6.53397e-14</value>
            <value variable="v" tolerance="1e-12">2.89187e-15</value>
	    <value variable="p" tolerance="2.3e-10">2.26904e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">5.31797e-14</value>
            <value variable="v" tolerance="1e-12">4.82522e-15</value>
	    <value variable="p" tolerance="1.1e-10">1.00154e-10</value>
        </metric>
    </metrics>
</test>

