<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>r-adaptation test</description>
    <executable>NekMesh</executable>
    <parameters>-v -m varopti:hyperelastic:scalingfile=2d_adapt.txt:maxiter=10:subiter=2:nq=2 2d_adapt.mcf 2d_adapt-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">2d_adapt.mcf</file>
        <file description="Input File 2">2d_adapt.geo</file>
        <file description="Input File 3">2d_adapt.txt</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Invalid at end: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
