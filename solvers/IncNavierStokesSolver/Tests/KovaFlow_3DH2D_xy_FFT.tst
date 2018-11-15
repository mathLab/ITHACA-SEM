<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 2D, flow in xy, FFT</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH2D_xy_FFT.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH2D_xy_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.13861e-08</value>
            <value variable="v" tolerance="1e-12">1.03921e-08</value>
            <value variable="w" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="1e-9">1.01794e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.88773e-07</value>
            <value variable="v" tolerance="1e-12">4.1301e-08</value>
            <value variable="w" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="1e-9">1.80295e-06</value>
        </metric>
    </metrics>
</test>
