<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=8, 16 Fourier modes, using implicit mapping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P8_16modes_Mapping-implicit.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P8_16modes_Mapping-implicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">1.46561e-05</value>
            <value variable="v" tolerance="1e-9">1.26223e-08</value>
            <value variable="w" tolerance="1e-9">1.08012e-07</value>
	    <value variable="p" tolerance="1e-9">1.08632e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.64559e-05</value>
            <value variable="v" tolerance="1e-9">2.00666e-07</value>
            <value variable="w" tolerance="1e-9">3.53111e-07</value>
	    <value variable="p" tolerance="1e-9">7.99001e-06</value>
        </metric>
    </metrics>
</test>
