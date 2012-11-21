<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Laminar Channel Flow 3D homogeneous 2D, P=3, 8x8 Fourier modes (MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH2D_MVM.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH2D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0</value>
            <value variable="v" tolerance="1e-6">6.55615e-16</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">1.23804e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">5.27402e-17</value>
            <value variable="v" tolerance="1e-6">1.05471e-15</value>
            <value variable="w" tolerance="1e-6">2.78544e-17</value>
            <value variable="p" tolerance="1e-6">2.99773e-15</value>
        </metric>
    </metrics>
</test>


