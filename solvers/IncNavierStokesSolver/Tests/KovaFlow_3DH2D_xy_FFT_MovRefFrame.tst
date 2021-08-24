<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 2D, flow in xy, MovingReferenceFrame, FFT</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH2D_xy_FFT_MovRefFrame.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH2D_xy_FFT_MovRefFrame.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="2.2e-5">0</value>
            <value variable="v" tolerance="2.23e-5">0</value>
            <value variable="w" tolerance="1e-9">0</value>
	    <value variable="p" tolerance="4.3e-5">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="7.1e-5">0</value>
            <value variable="v" tolerance="1.1e-4">0</value>
            <value variable="w" tolerance="1e-9">0</value>
            <value variable="p" tolerance="2.7e-4">0</value>
        </metric>
    </metrics>
</test>
