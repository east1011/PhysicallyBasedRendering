using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Audio;
using Microsoft.Xna.Framework.Content;
using Microsoft.Xna.Framework.GamerServices;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Input;
using Microsoft.Xna.Framework.Media;

namespace ILSBasicLightComplete
{

    public partial class Game1 : Microsoft.Xna.Framework.Game
    {
        GraphicsDeviceManager graphics;
        KeyboardState currentKeys;

        Effect effect;
        Texture2D texture;

        Matrix World, View, Projection;

        Vector3 lightDir = Vector3.Down;
        Vector3 lightPos = new Vector3(0, 15, 30);
        Vector3 spotPos = new Vector3(0, 30, 0);
        Vector3 spotDir = new Vector3(0, -1, -1);
        Vector3 cameraPos = new Vector3(0, 60, 60);

        public Game1()
        {
            graphics = new GraphicsDeviceManager(this);
            Content.RootDirectory = "Content";
        }

        protected override void Initialize()
        {
            graphics.PreferredBackBufferWidth = 1280;
            graphics.PreferredBackBufferHeight = 720;
            graphics.PreferMultiSampling = true;
            graphics.ApplyChanges();
            base.Initialize();
        }


        protected override void LoadContent()
        {
            effect = Content.Load<Effect>("DirectionalLight");

            texture = Content.Load<Texture2D>("grass");

            effect.Parameters["lightColor"].SetValue(Color.White.ToVector3());
            effect.Parameters["globalAmbient"].SetValue(Color.White.ToVector3());
            effect.Parameters["Ke"].SetValue(0.0f);
            effect.Parameters["Ka"].SetValue(0.2f);
            effect.Parameters["Kd"].SetValue(1.0f);
            effect.Parameters["Ks"].SetValue(0.3f);
            effect.Parameters["specularPower"].SetValue(100);
            effect.Parameters["spotPower"].SetValue(10);
            effect.Parameters["Texture"].SetValue(texture);


            World = Matrix.Identity;
            View = Matrix.CreateLookAt(cameraPos, Vector3.Zero, Vector3.Up);
            Projection = Matrix.CreatePerspectiveFieldOfView(MathHelper.PiOver4, GraphicsDevice.Viewport.AspectRatio, 1, 200);

            //Backreference calls
            CreateVertexBuffer();
            CreateIndexBuffer();
        }


        protected override void UnloadContent()
        {
        }


        protected override void Update(GameTime gameTime)
        {
            currentKeys = Keyboard.GetState();

            //Press Esc To Exit
            if (currentKeys.IsKeyDown(Keys.Escape))
                this.Exit();

            //Press Directional Keys to rotate camera
            if (currentKeys.IsKeyDown(Keys.Left))
            {
                cameraPos = Vector3.Transform(cameraPos, Matrix.CreateRotationY(-0.05f));
                View = Matrix.CreateLookAt(cameraPos, Vector3.Zero, Vector3.Up);
            }
            if (currentKeys.IsKeyDown(Keys.Right))
            {
                cameraPos = Vector3.Transform(cameraPos, Matrix.CreateRotationY(0.05f));
                View = Matrix.CreateLookAt(cameraPos, Vector3.Zero, Vector3.Up);
            }

            if(currentKeys.IsKeyDown(Keys.NumPad1)) effect.CurrentTechnique = effect.Techniques["DirectionalLight"];
            if(currentKeys.IsKeyDown(Keys.NumPad2)) effect.CurrentTechnique = effect.Techniques["PointLight"];
            if(currentKeys.IsKeyDown(Keys.NumPad3)) effect.CurrentTechnique = effect.Techniques["SpotLight"];

            lightPos = Vector3.Transform(lightPos, Matrix.CreateRotationY(0.03f));

            spotDir = Vector3.Transform(spotDir, Matrix.CreateRotationY(0.05f));
            spotDir.Normalize();

            lightDir = Vector3.Transform(lightDir,Matrix.CreateRotationZ(0.03f));
            lightDir.Normalize();


            base.Update(gameTime);
        }

        protected override void Draw(GameTime gameTime)
        {
            GraphicsDevice.Clear(Color.CornflowerBlue);

            GraphicsDevice.SetVertexBuffer(vertexBuffer);
            GraphicsDevice.Indices = indexBuffer;

            effect.Parameters["WVP"].SetValue(World * View * Projection);
            effect.Parameters["World"].SetValue(World);
            effect.Parameters["eyePosition"].SetValue(Vector3.Transform(cameraPos, Matrix.Invert(World)));
            if (effect.CurrentTechnique == effect.Techniques["SpotLight"])
            {
                effect.Parameters["lightDirection"].SetValue(Vector3.TransformNormal(spotDir, Matrix.Invert(World)));
                effect.Parameters["lightPosition"].SetValue(Vector3.Transform(spotPos, Matrix.Invert(World)));
            }
            else
            {
                effect.Parameters["lightDirection"].SetValue(Vector3.TransformNormal(lightDir, Matrix.Invert(World)));
                effect.Parameters["lightPosition"].SetValue(Vector3.Transform(lightPos, Matrix.Invert(World)));
            }
            effect.CurrentTechnique.Passes[0].Apply();

            GraphicsDevice.DrawIndexedPrimitives(PrimitiveType.TriangleList, 0, 0, number_of_vertices, 0, number_of_indices / 3);


            base.Draw(gameTime);
        }
    }
}
