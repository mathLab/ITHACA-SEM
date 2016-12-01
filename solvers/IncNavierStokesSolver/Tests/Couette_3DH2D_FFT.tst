<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Plane Couette Flow 3D homogeneous 2D, P=3, 4x4 Fourier modes (FFT)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I USEFFT=FFTW Couette_3DH2D.xml</parameters>
    <files>
        <file description="Session File">Couette_3DH2D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">0</value>
        </metric>
    </metrics>
</test>


